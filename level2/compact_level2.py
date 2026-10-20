#!/usr/bin/env python3
"""Create a review directory containing Level 2 summaries and best artifacts.

The predictor-specific raw trees can be very large.  This command copies the
best structure and the confidence/full-data files referenced by the result
manifest into flat, type-based directories before optionally moving raw run
directories into a reversible ``_archive`` directory.
"""

from __future__ import annotations

import argparse
import csv
import json
import shutil
import tarfile
from datetime import datetime, timezone
from pathlib import Path


def _copy(source: str, destination: Path) -> bool:
    path = Path(source) if source else None
    if not path or not path.is_file():
        return False
    destination.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(path, destination)
    return True


def _path_list(value: str) -> list[str]:
    """Read semicolon-separated paths emitted for Boltz full-data arrays."""

    return [item for item in str(value or "").split(";") if item]


def _model_state_dir(root: Path, model: str) -> Path:
    """Return a stable model/state directory for delivery artifacts."""
    if model.endswith("_apo"):
        return root / model[:-4] / "apo"
    if model.endswith("_holo"):
        return root / model[:-5] / "holo"
    return root / model / "single"


def _canonical_delivery_summary(run: Path, destination: Path) -> None:
    """Copy canonical metrics with paths rewritten into the delivery tree."""

    source = run / "metrics/canonical_metrics.csv"
    if not source.is_file():
        return
    rows: list[dict[str, str]] = []
    with source.open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle):
            design = row.get("design_id", "")
            model = row.get("model", "")
            structure = Path(row.get("structure_path", ""))
            if not structure.is_absolute():
                structure = (run / structure).resolve()
            if structure.is_file():
                target = _model_state_dir(destination / "structures", model) / f"{design}{structure.suffix or '.cif'}"
                row["structure_path"] = str(target.relative_to(destination)) if target.is_file() else ""
            else:
                row["structure_path"] = ""
            full_json = Path(row.get("full_data_path", ""))
            if not full_json.is_absolute():
                full_json = (run / full_json).resolve()
            json_target = _model_state_dir(destination / "confidences", model) / f"{design}_full_data.json"
            row["full_data_path"] = str(json_target.relative_to(destination)) if full_json.is_file() and json_target.is_file() else ""
            npz_targets: list[str] = []
            for item in _path_list(row.get("full_data_npz_path", "")):
                path = Path(item)
                if not path.is_absolute():
                    path = (run / path).resolve()
                target = _model_state_dir(destination / "confidences", model) / f"{design}_{path.name}"
                if path.is_file() and target.is_file():
                    npz_targets.append(str(target.relative_to(destination)))
            row["full_data_npz_path"] = ";".join(npz_targets)
            # The delivery summary itself is the authoritative metric table;
            # do not leave a stale absolute path to a run-internal summary JSON.
            row["summary_path"] = ""
            rows.append(row)
    if rows:
        destination.mkdir(parents=True, exist_ok=True)
        with (destination / "summary.tsv").open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
            writer.writeheader()
            writer.writerows(rows)
        metric_dir = destination / "metrics"
        metric_dir.mkdir(parents=True, exist_ok=True)
        with (metric_dir / "canonical_metrics.csv").open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)
        # Keep the JSON aggregate consistent with summary.tsv.  In
        # particular, do not copy run-internal absolute paths into review/.
        run_summary = run / "summary.json"
        if run_summary.is_file():
            try:
                payload = json.loads(run_summary.read_text(encoding="utf-8"))
                if isinstance(payload, dict):
                    payload["rows"] = rows
                    (destination / "summary.json").write_text(
                        json.dumps(payload, indent=2, ensure_ascii=False) + "\n",
                        encoding="utf-8",
                    )
            except (OSError, TypeError, ValueError):
                pass


def _validate_delivery(destination: Path) -> None:
    """Fail closed if the public handoff contains stale or absolute paths."""

    summary = destination / "summary.tsv"
    if not summary.is_file():
        raise RuntimeError("delivery summary.tsv is missing")
    required = {
        "design_id", "model", "status", "complex_iptm", "complex_ptm",
        "complex_mean_plddt", "complex_mean_pae", "complex_mean_interface_pae",
        "complex_ipsae", "chain_metrics_json", "chain_pair_metrics_json",
        "structure_path", "full_data_path", "full_data_npz_path", "usalign_rmsd",
        "usalign_chain_metrics_json",
    }
    with summary.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = set(reader.fieldnames or ())
        missing = sorted(required - fields)
        if missing:
            raise RuntimeError("delivery summary is missing columns: " + ", ".join(missing))
        rows = list(reader)
    if not rows:
        raise RuntimeError("delivery summary contains no rows")
    for row in rows:
        for field in ("structure_path", "full_data_path"):
            value = row.get(field, "")
            if value:
                path = Path(value)
                if path.is_absolute() or not (destination / path).is_file():
                    raise RuntimeError(f"invalid delivery path in {field}: {value}")
        for value in _path_list(row.get("full_data_npz_path", "")):
            path = Path(value)
            if path.is_absolute() or not (destination / path).is_file():
                raise RuntimeError(f"invalid delivery NPZ path: {value}")


def _pack_full_data(destination: Path) -> None:
    """Pack confidence/full-data artifacts and rewrite public references.

    The archive is the retained full-data artifact.  Individual JSON/NPZ files
    are removed only after the tarball has been closed and its members verified.
    """
    summary = destination / "summary.tsv"
    confidence_dir = destination / "confidences"
    if not summary.is_file() or not confidence_dir.is_dir():
        return
    files = sorted(path for path in confidence_dir.rglob("*") if path.is_file())
    if not files:
        return
    archive_dir = destination / "full_data"
    archive_dir.mkdir(parents=True, exist_ok=True)
    archive = archive_dir / "full_data.tar.gz"
    with tarfile.open(archive, "w:gz") as handle:
        for path in files:
            handle.add(path, arcname=str(path.relative_to(confidence_dir)))
    with tarfile.open(archive, "r:gz") as handle:
        members = {item.name for item in handle.getmembers() if item.isfile()}
    expected = {str(path.relative_to(confidence_dir)) for path in files}
    if members != expected:
        raise RuntimeError("full-data archive verification failed")

    rows: list[dict[str, str]] = []
    with summary.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = list(reader.fieldnames or [])
        if "full_data_archive_members" not in fields:
            fields.append("full_data_archive_members")
        for row in reader:
            members_for_row: list[str] = []
            for field in ("full_data_path", "full_data_npz_path"):
                raw = str(row.get(field, ""))
                if not raw:
                    continue
                for value in raw.split(";"):
                    rel = Path(value)
                    if rel.parts and rel.parts[0] == "confidences":
                        member = str(Path(*rel.parts[1:]))
                        if member in members:
                            members_for_row.append(member)
                if field == "full_data_npz_path" and raw:
                    row[field] = "full_data/full_data.tar.gz"
                elif field == "full_data_path" and raw:
                    row[field] = "full_data/full_data.tar.gz"
            row["full_data_archive_members"] = ";".join(sorted(set(members_for_row)))
            rows.append(row)
    with summary.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader(); writer.writerows(rows)
    metric = destination / "metrics/canonical_metrics.csv"
    if metric.is_file():
        with metric.open(encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle); metric_fields = list(reader.fieldnames or []); metric_rows = list(reader)
        if "full_data_archive_members" not in metric_fields:
            metric_fields.append("full_data_archive_members")
        by_key = {(r.get("design_id", ""), r.get("model", ""), r.get("state", "")): r for r in rows}
        for item in metric_rows:
            key = (item.get("design_id", ""), item.get("model", ""), item.get("state", ""))
            source = by_key.get(key, {})
            if source.get("full_data_path"):
                item["full_data_path"] = source["full_data_path"]
            if source.get("full_data_npz_path"):
                item["full_data_npz_path"] = source["full_data_npz_path"]
            item["full_data_archive_members"] = source.get("full_data_archive_members", "")
        with metric.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=metric_fields); writer.writeheader(); writer.writerows(metric_rows)
    shutil.rmtree(confidence_dir)


def compact(run: Path, destination: Path, *, archive: bool = False, delivery_only: bool = False) -> dict[str, int]:
    manifest = run / "results" / "result_manifest.tsv"
    if not manifest.is_file():
        raise FileNotFoundError(f"aggregate first: {manifest}")
    destination.mkdir(parents=True, exist_ok=True)
    structures = destination / "structures"
    confidences = destination / "confidences"
    copied_structures = copied_data = 0
    summary_rows: list[dict[str, str]] = []
    with manifest.open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            design_id = row.get("design_id", "").strip()
            if not design_id:
                continue
            for field, source in list(row.items()):
                if not field.endswith("_structure_path") or field in {"structure_path", "reference_structure_path", "backbone_structure_path"}:
                    continue
                model = field[: -len("_structure_path")]
                if not source:
                    continue
                suffix = Path(source).suffix or ".cif"
                target = _model_state_dir(structures, model) / f"{design_id}{suffix}"
                if _copy(source, target):
                    copied_structures += 1
                    row[field] = str(target.relative_to(destination))
                elif target.is_file():
                    row[field] = str(Path("structures") / target.name)
            for field, source in list(row.items()):
                if not field.endswith("_full_data_json"):
                    continue
                model = field[: -len("_full_data_json")]
                target = _model_state_dir(confidences, model) / f"{design_id}_full_data.json"
                if source and _copy(source, target):
                    copied_data += 1
                    row[field] = str(target.relative_to(destination))
                elif target.is_file():
                    row[field] = str(target.relative_to(destination))
            for field, source in list(row.items()):
                if not field.endswith("_output_dir") or "esmfold2" not in field:
                    continue
                model = field[: -len("_output_dir")]
                esmfold_dir = Path(source)
                if esmfold_dir.is_dir():
                    candidates = sorted(esmfold_dir.glob(f"{design_id}*_confidence.json"))
                    exact = esmfold_dir / f"{design_id}_confidence.json"
                    paths = [exact] if exact.is_file() else candidates[:1]
                    for path in paths:
                        target = _model_state_dir(confidences, model) / f"{design_id}_{path.name}"
                        if _copy(path.as_posix(), target):
                            copied_data += 1
            # The parser has already selected the top-ranked Boltz model and
            # records only its PAE/pLDDT NPZ paths.  Copy those exact paths;
            # scanning the raw directory would retain every diffusion sample.
            for field, value in list(row.items()):
                if not field.endswith("_full_data_npz"):
                    continue
                model = field[: -len("_full_data_npz")]
                for source in _path_list(value):
                    path = Path(source)
                    target = _model_state_dir(confidences, model) / f"{design_id}_{path.name}"
                    if _copy(path.as_posix(), target):
                        copied_data += 1
            summary_rows.append(row)

    if delivery_only and (run / "metrics/canonical_metrics.csv").is_file():
        # The delivery summary is metric-only.  The full result manifest stays
        # in the normal audit mode and is not copied into a delivery folder.
        _canonical_delivery_summary(run, destination)
    elif summary_rows:
        columns = list(summary_rows[0])
        with (destination / "summary.tsv").open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t")
            writer.writeheader()
            writer.writerows(summary_rows)

    if delivery_only and not (destination / "summary.tsv").is_file():
        # Never archive raw prediction artifacts when the metric handoff was
        # not produced successfully.
        raise RuntimeError("delivery-only closeout did not produce summary.tsv")

    copy_dirs = ("metrics",) if delivery_only else ("results", "metrics", "status")
    for name in copy_dirs:
        source = run / name
        if source.exists():
            shutil.copytree(source, destination / name, dirs_exist_ok=True)
    # copytree above restores the raw canonical table; rewrite it once more so
    # both public metric tables and summary.tsv use review-relative paths.
    if delivery_only and (run / "metrics/canonical_metrics.csv").is_file():
        _canonical_delivery_summary(run, destination)
    for name in ("summary.json", "ranking.csv"):
        source = run / name
        if source.is_file() and not (delivery_only and (destination / name).is_file()):
            shutil.copy2(source, destination / name)

    if delivery_only:
        _pack_full_data(destination)
        # Rebuild the JSON summary after packing so it carries archive paths.
        summary = destination / "summary.tsv"
        if summary.is_file() and (destination / "summary.json").is_file():
            try:
                payload = json.loads((destination / "summary.json").read_text(encoding="utf-8"))
                with summary.open(encoding="utf-8", newline="") as handle:
                    payload["rows"] = list(csv.DictReader(handle, delimiter="\t"))
                (destination / "summary.json").write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
            except (OSError, TypeError, ValueError):
                pass
        _validate_delivery(destination)

    if archive:
        stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
        archive_root = run / "_archive" / f"raw_{stamp}"
        archive_root.mkdir(parents=True, exist_ok=True)
        for name in ("artifacts", "inputs", "logs"):
            source = run / name
            if source.exists():
                shutil.move(str(source), str(archive_root / name))
    return {"structures": copied_structures, "confidence_files": copied_data}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", required=True, type=Path)
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument("--archive", action="store_true")
    parser.add_argument("--delivery-only", action="store_true", help="write only structures, confidence data, summary.tsv and metrics")
    args = parser.parse_args()
    counts = compact(args.run_dir.resolve(), args.out_dir.resolve(), archive=args.archive, delivery_only=args.delivery_only)
    print(f"compacted: {counts['structures']} structures, {counts['confidence_files']} confidence/full-data files")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
