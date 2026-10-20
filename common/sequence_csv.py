#!/usr/bin/env python3
"""Convert a sequence CSV into ProFold's normalized manifest contract.

The public entry point accepts CSV because that is convenient for design
libraries.  Internally we still use the strict manifest reader so paths are
resolved once and all downstream stages see the same row identity.
"""

from __future__ import annotations

import csv
from pathlib import Path
from typing import Mapping

from .manifest import normalize_row, write_tsv


def read_sequence_csv(path: str | Path) -> list[dict[str, str]]:
    source = Path(path).expanduser().resolve()
    if not source.is_file():
        raise FileNotFoundError(f"sequence CSV not found: {source}")
    with source.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        if not reader.fieldnames or "design_id" not in reader.fieldnames:
            raise ValueError("sequence CSV must contain design_id")
        if "sequence" not in reader.fieldnames and "fasta_path" not in reader.fieldnames:
            raise ValueError("sequence CSV must contain sequence or fasta_path")
        rows: list[dict[str, str]] = []
        for row in reader:
            if None in row or any(value is None for value in row.values()):
                raise ValueError(f"CSV line {reader.line_num}: column count differs from header")
            rows.append(normalize_row({str(k): str(v or "") for k, v in row.items()}, source))
    if not rows:
        raise ValueError(f"sequence CSV is empty: {source}")
    ids = [row["design_id"] for row in rows]
    duplicate = sorted({item for item in ids if ids.count(item) > 1})
    if duplicate:
        raise ValueError(f"duplicate design_id: {', '.join(duplicate)}")
    missing = [row["design_id"] for row in rows if not row.get("sequence") and not row.get("fasta_path")]
    if missing:
        raise ValueError(f"missing sequence/fasta_path for: {', '.join(missing[:10])}")
    return rows


def write_manifest_from_csv(csv_path: str | Path, manifest_path: str | Path) -> Path:
    rows = read_sequence_csv(csv_path)
    output = Path(manifest_path).expanduser().resolve()
    write_tsv(output, rows, preferred=("design_id", "sequence", "reference_structure_path"))
    return output
