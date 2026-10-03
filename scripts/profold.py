#!/usr/bin/env python3
"""Compact v38.2 facade for ProFold.

The legacy level1/level2 scripts remain the execution adapters.  This entry
point keeps the user contract small: one designs.tsv, one JSON config, and one
fixed summary.tsv per run.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from common.manifest import read_manifest, write_tsv  # noqa: E402
from common.pipeline import prepare_run  # noqa: E402
from report.generate_report import generate_report  # noqa: E402

SUMMARY_FIELDS = (
    "design_id", "stage", "status", "structure_path", "plddt", "ptm",
    "iptm", "interface_pae", "model_consensus", "selection", "error",
)


def _write_summary(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(f".{path.name}.tmp")
    with tmp.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(SUMMARY_FIELDS), delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows({field: row.get(field, "") for field in SUMMARY_FIELDS} for row in rows)
    tmp.replace(path)


def _load_json(path: str | Path) -> dict[str, object]:
    value = json.loads(Path(path).read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError(f"JSON object required: {path}")
    return value


def _canonical_rows(path: str | Path) -> list[dict[str, str]]:
    rows = read_manifest(path, require_sequence=True)
    source = Path(path).resolve().parent
    output: list[dict[str, str]] = []
    for row in rows:
        item = dict(row)
        mode = item.get("mode", "protein").strip().lower() or "protein"
        if mode not in {"protein", "cid"}:
            raise ValueError(f"{item['design_id']}: mode must be protein or cid")
        item["mode"] = mode
        # Compact aliases are deliberately translated once at the boundary.
        template = item.get("template", "").strip()
        if template and not item.get("target_spec_path"):
            item["target_spec_path"] = str((source / template).resolve()) if not Path(template).is_absolute() else template
        if mode == "cid":
            for key in ("cid_apo_json_path", "cid_holo_json_path"):
                if item.get(key) and not Path(item[key]).is_absolute():
                    item[key] = str((source / item[key]).resolve())
            item.setdefault("af3_mode", "cid")
        output.append(item)
    return output


def _fingerprint(input_path: Path, config_path: Path, level: str) -> str:
    h = hashlib.sha256()
    h.update(input_path.read_bytes())
    h.update(config_path.read_bytes())
    h.update(level.encode())
    h.update((ROOT / "VERSION").read_bytes())
    return h.hexdigest()


def prepare(args: argparse.Namespace) -> int:
    input_path = Path(args.input).resolve()
    config_path = Path(args.config).resolve()
    config = _load_json(config_path)
    rows = _canonical_rows(input_path)
    # write a normalized boundary manifest so downstream adapters never see
    # user-facing aliases or relative paths.
    normalized = Path(args.outdir).resolve() / "input.tsv"
    write_tsv(normalized, rows, preferred=("design_id", "sequence", "mode", "template", "ligand", "reference_structure"))
    level = args.level.lower()
    if level not in {"level1", "level2"}:
        raise ValueError("--level must be level1 or level2")
    if any(row.get("mode") == "cid" for row in rows):
        config = dict(config)
        config["af3_mode"] = "cid"
    out, _ = prepare_run(normalized, args.outdir, level=level, shard_size=args.shard_size, config=config)
    fp = _fingerprint(out / "input.tsv", out / "config.json", level)
    (out / "fingerprint.json").write_text(json.dumps({"sha256": fp, "level": level}, indent=2) + "\n", encoding="utf-8")
    print(f"prepared {level}: {len(rows)} designs -> {out}")
    return 0


def _run_script(level: str, action: str) -> Path:
    return ROOT / level / f"{action}_{level}.py"


def run(args: argparse.Namespace) -> int:
    run_dir = Path(args.run_dir).resolve()
    config = json.loads((run_dir / "config.json").read_text(encoding="utf-8"))
    fingerprint_path = run_dir / "fingerprint.json"
    if fingerprint_path.is_file():
        stored = _load_json(fingerprint_path).get("sha256", "")
        current = _fingerprint(run_dir / "input.tsv", run_dir / "config.json", args.level)
        if stored != current:
            raise SystemExit("resume refused: input/config/version fingerprint changed")
    n = int((run_dir / "tasks/n_shards.txt").read_text().strip())
    worker = ROOT / ("level1" if args.level == "level1" else "level2") / ("run_level1.py" if args.level == "level1" else "run_level2.py")
    for index in range(n):
        command = [sys.executable, str(worker), "--run-dir", str(run_dir), "--shard-index", str(index)]
        if args.no_resume:
            command.append("--no-resume")
        subprocess.run(command, check=True)
    return 0


def _summary_rows(run_dir: Path, level: str) -> list[dict[str, str]]:
    from common.manifest import read_manifest
    rows = read_manifest(run_dir / "results/result_manifest.tsv")
    models = ("esmfold", "af3") if level == "level1" else ("protenix", "boltz2", "opendde")
    stateful = level == "level2" and any(str(row.get("states", "")).strip() for row in rows)
    metric_models = tuple(f"{m}_holo" for m in models) if stateful else models
    out = []
    for row in rows:
        success = [m for m in metric_models if row.get(f"{m}_status") == "success"]
        preferred = ("af3", "esmfold") if level == "level1" else metric_models
        ordered = [m for m in preferred if m in success]
        structure = next((row.get(f"{m}_structure_path", "") for m in ordered if row.get(f"{m}_structure_path")), "")
        plddt = next((row.get(f"{m}_mean_plddt", "") for m in ordered if row.get(f"{m}_mean_plddt", "") != ""), "")
        ptm = next((row.get(f"{m}_ptm", "") for m in ordered if row.get(f"{m}_ptm", "") != ""), "")
        iptm = next((row.get(f"{m}_iptm", "") for m in ordered if row.get(f"{m}_iptm", "") != ""), "")
        out.append({
            "design_id": row.get("design_id", ""), "stage": level,
            "status": row.get("status", "missing"), "structure_path": structure,
            "plddt": plddt, "ptm": ptm, "iptm": iptm, "interface_pae": "",
            "model_consensus": f"{len(success)}/{len(metric_models)}",
            "selection": "pass" if row.get("status") == "success" else "",
            "error": row.get("error", ""),
        })
    return out


def aggregate(args: argparse.Namespace) -> int:
    run_dir = Path(args.run_dir).resolve()
    level_dir = ROOT / args.level
    script = level_dir / f"aggregate_{args.level}.py"
    subprocess.run([sys.executable, str(script), "--run-dir", str(run_dir)], check=True)
    rows = _summary_rows(run_dir, args.level)
    _write_summary(run_dir / "summary.tsv", rows)
    print(f"summary: {run_dir / 'summary.tsv'} ({len(rows)} designs)")
    return 0


def report(args: argparse.Namespace) -> int:
    output = generate_report(
        args.run_dir,
        top_n=args.top_n,
        render_structures=not args.no_structures,
        renderer=args.renderer,
        resolution=args.render_resolution,
        aa_samples=args.aa_samples,
        pdf_name=args.pdf_name,
    )
    print(f"report: {output}")
    return 0


def main() -> int:
    parser = argparse.ArgumentParser(description="ProFold v38.2 compact facade")
    sub = parser.add_subparsers(dest="command", required=True)
    p = sub.add_parser("prepare"); p.add_argument("--input", required=True); p.add_argument("--config", required=True); p.add_argument("--outdir", required=True); p.add_argument("--level", default="level1"); p.add_argument("--shard-size", type=int, default=1); p.set_defaults(func=prepare)
    p = sub.add_parser("run"); p.add_argument("--run-dir", required=True); p.add_argument("--level", required=True, choices=("level1", "level2")); p.add_argument("--no-resume", action="store_true"); p.set_defaults(func=run)
    p = sub.add_parser("aggregate"); p.add_argument("--run-dir", required=True); p.add_argument("--level", required=True, choices=("level1", "level2")); p.set_defaults(func=aggregate)
    p = sub.add_parser("report", help="generate plots, Top-N structure panels, and a PDF report")
    p.add_argument("--run-dir", required=True)
    p.add_argument("--top-n", type=int, default=10)
    p.add_argument("--no-structures", action="store_true")
    p.add_argument("--renderer", default="/xcfhome/yhliu/14_magpcr/z_zcodex/pipeline_scripts/visualization/render_structure_panel.py")
    p.add_argument("--render-resolution", default="1200x900")
    p.add_argument("--aa-samples", type=int, default=4)
    p.add_argument("--pdf-name", default="profold_level2_report.pdf")
    p.set_defaults(func=report)
    args = parser.parse_args()
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
