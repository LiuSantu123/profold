#!/usr/bin/env python3
"""Batch interface for optional Rosetta/BoltzGen-style structure filters.

Rosetta installations differ by site, so ProFold does not embed a fragile
command line.  A site supplies one executable/script that accepts the input
TSV and output JSON path.  The controller invokes it once per run, avoiding
one Rosetta startup per design.
"""

from __future__ import annotations

import json
import math
import subprocess
import csv
from pathlib import Path


FILTER_FIELDS = (
    "rosetta_interface_complementarity",
    "rosetta_sasa",
    "rosetta_delta_sasa",
    "rosetta_radius_of_gyration",
)


def _validate_payload(payload: object, design_ids: set[str]) -> list[dict[str, object]]:
    """Validate and normalize the batch script's JSON contract."""

    if isinstance(payload, dict):
        rows = payload.get("rows")
    else:
        rows = payload
    if not isinstance(rows, list):
        raise ValueError("Rosetta output must be a JSON list or an object with a 'rows' list")
    normalized: list[dict[str, object]] = []
    for item in rows:
        if not isinstance(item, dict) or not str(item.get("design_id", "")):
            raise ValueError("each Rosetta row must contain design_id")
        row = {str(key): value for key, value in item.items()}
        for field in FILTER_FIELDS:
            if field not in row or row[field] in (None, ""):
                continue
            try:
                value = float(row[field])
            except (TypeError, ValueError) as exc:
                raise ValueError(f"{field} must be numeric for {row['design_id']}") from exc
            if not math.isfinite(value):
                raise ValueError(f"{field} must be finite for {row['design_id']}")
            row[field] = value
        normalized.append(row)
    unknown = sorted({str(row["design_id"]) for row in normalized} - design_ids)
    if unknown:
        raise ValueError("Rosetta output contains unknown design_id(s): " + ", ".join(unknown))
    return normalized


def run_rosetta_filters(run_dir: str | Path, script: str | Path, *, timeout: int = 3600) -> Path:
    run = Path(run_dir).resolve()
    manifest = run / "results" / "result_manifest.tsv"
    output = run / "metrics" / "rosetta_filters.json"
    output.parent.mkdir(parents=True, exist_ok=True)
    command = [str(Path(script).expanduser()), "--input", str(manifest), "--output", str(output)]
    try:
        completed = subprocess.run(command, cwd=str(run), capture_output=True, text=True, timeout=timeout)
        if completed.returncode != 0:
            payload = {"status": "failed", "error": (completed.stderr or completed.stdout)[-2000:]}
            output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
        elif not output.is_file():
            output.write_text(json.dumps({"status": "failed", "error": "Rosetta script did not create output JSON"}, indent=2) + "\n", encoding="utf-8")
        else:
            with output.open(encoding="utf-8") as handle:
                payload = json.load(handle)
            with manifest.open(encoding="utf-8", newline="") as handle:
                design_ids = {str(row.get("design_id", "")) for row in csv.DictReader(handle, delimiter="\t")}
            normalized = _validate_payload(payload, design_ids)
            output.write_text(json.dumps(normalized, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    except Exception as exc:
        output.write_text(json.dumps({"status": "failed", "error": f"{type(exc).__name__}: {exc}"}, indent=2) + "\n", encoding="utf-8")
    return output
