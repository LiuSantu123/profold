#!/usr/bin/env python3
"""Run Rosetta InterfaceAnalyzer on the AF3 structure for each design.

The script is intentionally a batch adapter for ``rosetta_metrics.py``.  It
accepts a result manifest and emits one JSON row per design, using AF3's
structure when available and the predictor best structure otherwise.
"""
from __future__ import annotations

import argparse
import csv
import json
import re
import subprocess
import tempfile
from pathlib import Path


def _chains(path: Path) -> list[str]:
    try:
        import gemmi
        structure = gemmi.read_structure(str(path))
        return [str(chain.name) for chain in structure[0]]
    except Exception:
        return ["A", "B"]


def _score(path: Path) -> dict[str, float]:
    lines = path.read_text(errors="replace").splitlines()
    headers: list[str] = []
    values: list[str] = []
    for line in lines:
        if line.startswith("SEQUENCE:"):
            continue
        fields = line.split()
        if "dSASA_int" in fields and "total_score" in fields:
            headers = fields
        elif headers and len(fields) >= len(headers) and fields[0] not in {"SCORE:", "description"}:
            values = fields[:len(headers)]
    if not headers or not values:
        return {}
    raw = dict(zip(headers, values))
    aliases = {
        "rosetta_interface_complementarity": ("sc_value", "shape_complementarity"),
        "rosetta_sasa": ("complexed_sasa",),
        "rosetta_delta_sasa": ("dSASA_int",),
        "rosetta_radius_of_gyration": ("rg",),
    }
    result: dict[str, float] = {}
    for target, names in aliases.items():
        for name in names:
            try:
                result[target] = float(raw[name])
                break
            except (KeyError, TypeError, ValueError):
                continue
    return result


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True, type=Path)
    ap.add_argument("--output", required=True, type=Path)
    ap.add_argument("--executable", required=True, type=Path)
    ap.add_argument("--timeout", type=int, default=600)
    args = ap.parse_args()
    rows: dict[str, dict[str, object]] = {}
    with args.input.open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            design = row.get("design_id", "")
            if design in rows:
                continue
            path = Path(row.get("af3_structure_path", "") or "")
            if not path.is_file():
                for key, value in row.items():
                    if key.endswith("_structure_path") and value and Path(value).is_file():
                        path = Path(value); break
            item: dict[str, object] = {"design_id": design, "structure_path": str(path)}
            if not path.is_file():
                item["status"] = "missing"; rows[design] = item; continue
            chain_ids = _chains(path)
            interface = "_".join(chain_ids[:2]) if len(chain_ids) >= 2 else "A_B"
            with tempfile.TemporaryDirectory(prefix="profold_rosetta_") as tmp:
                score = Path(tmp) / "score.sc"
                cmd = [str(args.executable), "-s", str(path), "-interface", interface,
                       "-pack_input", "-out:file:scorefile", str(score), "-overwrite"]
                try:
                    proc = subprocess.run(cmd, capture_output=True, text=True,
                                          timeout=args.timeout, cwd=tmp)
                    if proc.returncode:
                        item.update(status="failed", error=(proc.stderr or proc.stdout)[-1000:])
                    else:
                        item.update(status="success", **_score(score))
                except Exception as exc:
                    item.update(status="failed", error=f"{type(exc).__name__}: {exc}")
            rows[design] = item
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(list(rows.values()), indent=2) + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
