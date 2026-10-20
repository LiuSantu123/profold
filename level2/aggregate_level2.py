#!/usr/bin/env python3
"""Aggregate Level 2 shard outputs."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

V37 = Path(__file__).resolve().parents[1]
if str(V37) not in sys.path:
    sys.path.insert(0, str(V37))

from common.manifest import write_csv, write_tsv  # noqa: E402
from common.pipeline import join_results, metric_rows, status_rows  # noqa: E402
from common.pipeline import load_config  # noqa: E402
from level2.run_level2 import configured_models  # noqa: E402
from level2.align_level2 import align_run, DEFAULT_USALIGN  # noqa: E402
from common.metric_schema import CANONICAL_FIELDS, canonical_metric_row  # noqa: E402
from common.runner import atomic_json  # noqa: E402
from level2.rosetta_metrics import run_rosetta_filters  # noqa: E402
from level2.structure_filters import run_geometry_filters  # noqa: E402
import json  # noqa: E402
import csv  # noqa: E402
from tools.monitors.esmfold2_monitor import scan as scan_esmfold2, write_csv as write_esmfold2_csv  # noqa: E402


def _write_stage_timings(run: Path) -> None:
    """Collect per-shard predictor timings into one tuning-friendly table."""

    rows: list[dict[str, object]] = []
    for path in sorted((run / "artifacts").glob("shard_*/stage_timings.json")):
        try:
            payload = json.loads(path.read_text(encoding="utf-8"))
        except (OSError, TypeError, ValueError):
            continue
        if isinstance(payload, list):
            rows.extend(item for item in payload if isinstance(item, dict))
    if rows:
        write_tsv(
            run / "metrics/stage_timings.tsv",
            rows,
            preferred=(
                "shard_index", "model", "state", "design_count", "result_count",
                "success_count", "elapsed_seconds", "resumed",
            ),
            include_core=False,
        )


def _cid_enabled(config: dict[str, object]) -> bool:
    mode = str(config.get("level2_mode", config.get("cid_mode", "single"))).strip().lower()
    return mode in {"cid", "apo_holo", "holo_apo"}


def _number(value: object) -> float | None:
    try:
        number = float(value)
        return number if number == number else None
    except (TypeError, ValueError):
        return None


def _cid_delta_rows(rows: list[dict[str, str]], models: tuple[str, ...]) -> list[dict[str, object]]:
    metrics = (
        "iptm", "ptm", "ranking_score", "mean_plddt", "overall_ipae",
        "overall_ipsae_mean", "overall_avg_if_pae", "confidence_score",
        "complex_iplddt", "affinity_pred_value", "affinity_probability_binary",
    )
    output: list[dict[str, object]] = []
    for row in rows:
        for model in models:
            item: dict[str, object] = {
                "design_id": row.get("design_id", ""),
                "model": model,
                "apo_status": row.get(f"{model}_apo_status", "missing"),
                "holo_status": row.get(f"{model}_holo_status", "missing"),
            }
            for metric in metrics:
                apo = _number(row.get(f"{model}_apo_{metric}"))
                holo = _number(row.get(f"{model}_holo_{metric}"))
                item[f"apo_{metric}"] = "" if apo is None else apo
                item[f"holo_{metric}"] = "" if holo is None else holo
                item[f"delta_{metric}"] = "" if apo is None or holo is None else holo - apo
            output.append(item)
    return output


def _merge_post_metrics(canonical_rows: list[dict[str, object]], run: Path) -> None:
    """Attach CPU post-processing values without changing predictor rows."""
    ref = run / "metrics/usalign_reference.csv"
    ref_rows: dict[tuple[str, str], dict[str, str]] = {}
    if ref.is_file():
        with ref.open(encoding="utf-8", newline="") as handle:
            for item in csv.DictReader(handle):
                key = (str(item.get("design_id", "")), str(item.get("model", "")))
                ref_rows[key] = dict(item)
        for row in canonical_rows:
            item = ref_rows.get((str(row.get("design_id", "")), str(row.get("model", ""))))
            if item:
                row["usalign_rmsd"] = item.get("rmsd", "")
                row["usalign_tm_score"] = item.get("tm_score_norm_structure1", "")
    chain_file = run / "metrics/usalign_detailed/chain_pairs.csv"
    if chain_file.is_file():
        grouped: dict[str, list[dict[str, object]]] = {}
        with chain_file.open(encoding="utf-8", newline="") as handle:
            for item in csv.DictReader(handle):
                design = str(item.get("design_id", ""))
                grouped.setdefault(design, []).append({
                    "backbone_chain": item.get("backbone_chain", ""),
                    "predicted_chain": item.get("predicted_chain", ""),
                    "status": item.get("status", ""),
                    "rmsd": item.get("rmsd", ""),
                    "tm_score_norm_predicted": item.get("tm_score_norm_predicted", ""),
                    "tm_score_norm_backbone": item.get("tm_score_norm_backbone", ""),
                    "aligned_length": item.get("aligned_length", ""),
                })
        for row in canonical_rows:
            detail_id = f"{row.get('design_id', '')}__model={row.get('model', '')}"
            pairs = grouped.get(detail_id, [])
            if pairs:
                row["usalign_chain_metrics_json"] = json.dumps(pairs, ensure_ascii=False, sort_keys=True)
    rosetta = run / "metrics/rosetta_filters.json"
    if rosetta.is_file():
        try:
            payload = json.loads(rosetta.read_text(encoding="utf-8"))
            items = payload if isinstance(payload, list) else payload.get("rows", [])
            by_id = {str(item.get("design_id")): item for item in items if isinstance(item, dict)}
            for row in canonical_rows:
                values = by_id.get(str(row.get("design_id")), {})
                for field in ("rosetta_interface_complementarity", "rosetta_sasa", "rosetta_delta_sasa", "rosetta_radius_of_gyration", "rosetta_isoelectric_point"):
                    if field in values:
                        row[field] = values[field]
        except (OSError, TypeError, ValueError):
            pass
    geometry = run / "metrics/geometry_filters.json"
    if geometry.is_file():
        try:
            items = json.loads(geometry.read_text(encoding="utf-8"))
            for item in items if isinstance(items, list) else []:
                for row in canonical_rows:
                    same_design = row.get("design_id") == item.get("design_id")
                    same_model = row.get("model") == item.get("model")
                    # AF3 geometry is the preferred physical-filter source;
                    # propagate it to each Level2 predictor row for the same
                    # design so screening is model-comparable.
                    af3_item = item.get("model") == "af3"
                    if same_design and (same_model or af3_item):
                        row["geometry_radius_of_gyration"] = item.get("radius_of_gyration", "")
                        row["geometry_sasa"] = item.get("sasa", "")
                        row["geometry_delta_sasa"] = item.get("delta_sasa", "")
                        row["geometry_isoelectric_point"] = item.get("isoelectric_point", "")
                        row["geometry_source_structure"] = item.get("structure_path", "")
                        row["geometry_isoelectric_point"] = item.get("isoelectric_point", "")
                        row["geometry_interface_pairs_json"] = json.dumps(item.get("interface_pairs", {}), ensure_ascii=False, sort_keys=True)
                        pairs = item.get("interface_pairs", {})
                        if isinstance(pairs, dict):
                            densities = [value.get("contact_density") for value in pairs.values() if isinstance(value, dict) and isinstance(value.get("contact_density"), (int, float))]
                            if densities:
                                row["geometry_interface_contact_density"] = sum(densities) / len(densities)
        except (OSError, TypeError, ValueError):
            pass


def main() -> int:
    parser = argparse.ArgumentParser(description="Aggregate ProFold Level 2")
    parser.add_argument("--run-dir", required=True)
    args = parser.parse_args()
    run = Path(args.run_dir).resolve()
    base_models = configured_models(load_config(run))
    _write_stage_timings(run)
    cid = _cid_enabled(load_config(run))
    models = tuple(
        f"{model}_{state}" for model in base_models for state in ("apo", "holo")
    ) if cid else base_models
    rows = join_results(run, model_names=models)
    write_tsv(run / "results/result_manifest.tsv", rows)
    write_csv(run / "metrics/level2.csv", metric_rows(rows, models))
    canonical_rows = [
        canonical_metric_row(row, model=model, state=(model.rsplit("_", 1)[1] if cid else "single"))
        for row in rows for model in models
    ]
    # Written again after USalign/Rosetta so the public table contains all
    # post-prediction metrics while retaining the same stable columns.
    write_csv(run / "metrics/canonical_metrics.csv", canonical_rows, preferred=CANONICAL_FIELDS, include_core=False)
    write_tsv(run / "status/status.tsv", status_rows(rows, "level2", models))
    if cid:
        write_csv(run / "metrics/level2_cid_delta.csv", _cid_delta_rows(rows, base_models))
        for state in ("apo", "holo"):
            # Restrict the monitor to the state directory; this prevents
            # duplicate design IDs in the two-state reports.
            esmfold2_rows = scan_esmfold2(run / "artifacts")
            state_rows = [
                row for row in esmfold2_rows
                if f"/models/esmfold2/{state}/" in str(row.get("structure_path", ""))
            ]
            if state_rows:
                write_esmfold2_csv(state_rows, run / "metrics" / f"esmfold2_monitor_{state}.csv")
    else:
        esmfold2_rows = scan_esmfold2(run / "artifacts")
        if esmfold2_rows:
            write_esmfold2_csv(esmfold2_rows, run / "metrics/esmfold2_monitor.csv")
    # Structural comparison is a post-prediction stage.  It always writes the
    # internal-model table; reference rows are populated when the input
    # manifest records reference_structure_path/backbone_structure_path.
    usalign_path = str(load_config(run).get("usalign_path", DEFAULT_USALIGN))
    try:
        alignment = align_run(
            run,
            usalign_path=usalign_path,
            workers=int(load_config(run).get("usalign_workers", 4)),
        )
        print("USalign: " + ", ".join(f"{name}={path}" for name, path in alignment.items()))
    except Exception as exc:
        # Do not hide prediction results if an optional alignment executable is
        # unavailable; the alignment tables contain explicit missing status.
        print(f"USalign alignment skipped: {type(exc).__name__}: {exc}")
    rosetta_script = load_config(run).get("rosetta_filter_script")
    geometry_output = run_geometry_filters(run)
    print(f"Geometry filters: {geometry_output}")
    if rosetta_script:
        output = run_rosetta_filters(
            run,
            str(rosetta_script),
            timeout=int(load_config(run).get("rosetta_filter_timeout", 3600)),
        )
        print(f"Rosetta filters: {output}")
    _merge_post_metrics(canonical_rows, run)
    write_csv(run / "metrics/canonical_metrics.csv", canonical_rows, preferred=CANONICAL_FIELDS, include_core=False)
    counts: dict[str, int] = {}
    for row in rows:
        counts[str(row.get("status", "missing"))] = counts.get(str(row.get("status", "missing")), 0) + 1
    # Run-level metric aggregates are deliberately generated before closeout
    # moves raw artifacts.  They are small, model-neutral, and sufficient for
    # ranking without reopening predictor-specific output trees.
    atomic_json(
        run / "summary.json",
        {"rows": canonical_rows, "status_counts": counts, "models": list(models)},
    )
    ranking_rows = []
    for row in canonical_rows:
        score = _number(row.get("complex_ipsae"))
        if score is None:
            score = _number(row.get("complex_iptm"))
        if score is None:
            score = _number(row.get("complex_ptm"))
        ranking_rows.append({
            "design_id": row.get("design_id", ""),
            "model": row.get("model", ""),
            "state": row.get("state", "single"),
            "status": row.get("status", ""),
            "ranking_score": "" if score is None else score,
            "complex_ipsae": row.get("complex_ipsae", ""),
            "complex_iptm": row.get("complex_iptm", ""),
            "complex_ptm": row.get("complex_ptm", ""),
            "complex_mean_plddt": row.get("complex_mean_plddt", ""),
        })
    ranking_rows.sort(key=lambda item: (_number(item.get("ranking_score")) is not None, _number(item.get("ranking_score")) or float("-inf")), reverse=True)
    write_csv(
        run / "ranking.csv",
        ranking_rows,
        preferred=("design_id", "model", "state", "status", "ranking_score"),
        include_core=False,
    )
    print(f"aggregated level2: {len(rows)} rows -> {run / 'results/result_manifest.tsv'}; {counts}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
