#!/usr/bin/env python3
"""Run USalign on Level 2 structures.

The alignment stage is deliberately post-prediction: it reads the aggregated
manifest, compares every pair of available model structures for a design, and
compares each model to the design's reference structure when one is recorded.
"""
from __future__ import annotations

import csv
import itertools
import json
import subprocess
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from typing import Iterable

DEFAULT_USALIGN = "/data/backup/softwares/alphafold3_v4/03_tools/USalign/USalign"


def _write(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields: list[str] = []
    for row in rows:
        for key in row:
            if key not in fields:
                fields.append(key)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields or ["status"], extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def _parse(text: str, returncode: int) -> dict[str, object]:
    import re
    tm = [float(x) for x in re.findall(r"TM-score=\s*([0-9.eE+-]+)", text)]
    rmsd = re.search(r"Aligned length=\s*(\d+),\s*RMSD=\s*([0-9.eE+-]+)", text)
    ident = re.search(r"Seq_ID=n_identical/n_aligned=\s*([0-9.eE+-]+)", text)
    return {
        "status": "success" if returncode == 0 and len(tm) >= 2 else "failed",
        "aligned_length": int(rmsd.group(1)) if rmsd else "",
        "rmsd": float(rmsd.group(2)) if rmsd else "",
        "seq_identity": float(ident.group(1)) if ident else "",
        "tm_score_norm_structure1": tm[0] if len(tm) > 0 else "",
        "tm_score_norm_structure2": tm[1] if len(tm) > 1 else "",
        "error": "" if returncode == 0 and len(tm) >= 2 else f"USalign returncode={returncode}",
    }


def _run(a: Path, b: Path, executable: Path, timeout: int) -> dict[str, object]:
    if not a.is_file() or not b.is_file():
        return {"status": "missing", "error": f"missing structure: {a if not a.is_file() else b}"}
    if not executable.is_file():
        return {"status": "missing", "error": f"USalign not found: {executable}"}
    try:
        proc = subprocess.run([str(executable), str(a), str(b), "-mol", "prot", "-ter", "2"],
                              capture_output=True, text=True, timeout=timeout)
    except Exception as exc:
        return {"status": "failed", "error": f"{type(exc).__name__}: {exc}"}
    return _parse(proc.stdout or proc.stderr, proc.returncode)


def _reference(row: dict[str, str], base: Path) -> Path | None:
    for key in ("reference_structure_path", "backbone_structure_path", "reference_path", "backbone_path"):
        value = row.get(key, "").strip()
        if value:
            path = Path(value).expanduser()
            return path if path.is_absolute() else (base / path).resolve()
    return None


def align_run(run_dir: str | Path, *, usalign_path: str | Path = DEFAULT_USALIGN,
              timeout: int = 300, keep_raw: bool = False, workers: int = 1) -> dict[str, Path]:
    run = Path(run_dir).resolve()
    manifest = run / "results" / "result_manifest.tsv"
    if not manifest.is_file():
        raise FileNotFoundError(f"aggregate first: {manifest}")
    with manifest.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    exe = Path(usalign_path).expanduser().resolve()
    internal: list[dict[str, object]] = []
    reference: list[dict[str, object]] = []
    internal_jobs: list[tuple[dict[str, object], Path, Path]] = []
    reference_jobs: list[tuple[dict[str, object], Path, Path]] = []
    for row in rows:
        design = row.get("design_id", "")
        models: list[tuple[str, Path]] = []
        for key, value in row.items():
            if not key.endswith("_structure_path") or not value:
                continue
            model = key[:-len("_structure_path")]
            if model in {"structure", "reference", "backbone"}:
                continue
            path = Path(value)
            if not path.is_absolute():
                path = (run / path).resolve()
            if path.is_file():
                models.append((model, path))
        for (m1, p1), (m2, p2) in itertools.combinations(models, 2):
            internal_jobs.append(({"design_id": design, "model1": m1, "model2": m2,
                                  "structure1_path": str(p1), "structure2_path": str(p2)}, p1, p2))
        ref = _reference(row, run)
        if ref:
            for model, pred in models:
                reference_jobs.append(({"design_id": design, "model": model,
                                        "reference_structure_path": str(ref),
                                        "predicted_structure_path": str(pred)}, pred, ref))
    def run_job(job: tuple[dict[str, object], Path, Path]) -> dict[str, object]:
        metadata, left, right = job
        return {**metadata, **_run(left, right, exe, timeout)}
    if workers < 1:
        raise ValueError("workers must be >= 1")
    if workers == 1:
        internal = [run_job(job) for job in internal_jobs]
        reference = [run_job(job) for job in reference_jobs]
    else:
        with ThreadPoolExecutor(max_workers=workers) as pool:
            internal = list(pool.map(run_job, internal_jobs))
            reference = list(pool.map(run_job, reference_jobs))
    out_internal = run / "metrics" / "usalign_internal.csv"
    out_reference = run / "metrics" / "usalign_reference.csv"
    _write(out_internal, internal)
    _write(out_reference, reference)
    detailed_dir = run / "metrics" / "usalign_detailed"
    detailed_summary = ""
    if reference:
        # Reuse the existing chain-aware implementation for reference pairs;
        # the lightweight table above remains fast for model-model screening.
        detailed_input = run / "metrics" / "usalign_reference_input.csv"
        with detailed_input.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=("design_id", "reference_structure_path", "predicted_structure_path"), extrasaction="ignore")
            writer.writeheader()
            # structure_align intentionally rejects duplicate design IDs.  A
            # design can have several model-vs-reference rows, so make the
            # temporary IDs unique while retaining the original public table.
            for item in reference:
                detail = dict(item)
                detail["design_id"] = f"{item.get('design_id', '')}__model={item.get('model', '')}"
                writer.writerow(detail)
        try:
            from structure_align.align_structures import align_rows
            paths = align_rows(
                detailed_input,
                detailed_dir,
                usalign_path=exe,
                workers=max(1, workers),
                timeout=timeout,
                keep_raw=keep_raw,
            )
            detailed_summary = str(paths["summary"])
        except Exception as exc:
            (detailed_dir / "error.txt").parent.mkdir(parents=True, exist_ok=True)
            (detailed_dir / "error.txt").write_text(f"{type(exc).__name__}: {exc}\n", encoding="utf-8")
    summary = run / "metrics" / "usalign_summary.json"
    summary.write_text(json.dumps({"internal_pairs": len(internal), "reference_pairs": len(reference),
                                   "usalign": str(exe), "detailed_summary": detailed_summary}, indent=2) + "\n", encoding="utf-8")
    return {"internal": out_internal, "reference": out_reference, "summary": summary}


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-dir", required=True)
    parser.add_argument("--usalign", default=DEFAULT_USALIGN)
    parser.add_argument("--timeout", type=int, default=300)
    parser.add_argument("--workers", type=int, default=1)
    args = parser.parse_args()
    for name, path in align_run(args.run_dir, usalign_path=args.usalign, timeout=args.timeout, workers=args.workers).items():
        print(f"{name}: {path}")
