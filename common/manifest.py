#!/usr/bin/env python3
"""Stable, audit-friendly TSV manifest helpers.

The manifest is the only interface between levels.  Model-specific columns are
allowed, but ``design_id`` is always the stable primary key.
"""

from __future__ import annotations

import csv
import hashlib
import json
import re
from pathlib import Path
from typing import Iterable, Mapping


DESIGN_ID_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.-]*$")
PATH_FIELDS = {
    "fasta_path",
    "structure_path",
    "esmfold_structure_path",
    "af3_json_path",
    "cid_apo_json_path",
    "cid_holo_json_path",
    "target_spec_path",
    "protenix_input_path",
    "boltz_input_path",
    "opendde_input_path",
    "model_input_path",
    "source_manifest",
    "af3_structure_path",
    "protenix_structure_path",
    "boltz2_structure_path",
    "opendde_structure_path",
}
CORE_FIELDS = [
    "design_id",
    "parent_design_id",
    "sequence",
    "fasta_path",
    "structure_path",
    "source_manifest",
    "target_spec_path",
    "af3_json_path",
    "protenix_input_path",
    "boltz_input_path",
    "model_input_path",
    "input_kind",
    "target_id",
    "ligand_id",
    "chain_map",
    "metadata_json",
]


def _string(value: object) -> str:
    if value is None:
        return ""
    if isinstance(value, (dict, list, tuple)):
        return json.dumps(value, ensure_ascii=False, sort_keys=True)
    return str(value)


def _resolve(value: str, base_dir: Path) -> str:
    if not value:
        return ""
    path = Path(value).expanduser()
    if not path.is_absolute():
        path = (base_dir / path).resolve()
    return str(path)


def read_fasta(path: str | Path) -> list[tuple[str, str]]:
    """Read a small FASTA without importing optional bioinformatics packages."""

    records: list[tuple[str, str]] = []
    name: str | None = None
    chunks: list[str] = []
    for raw in Path(path).read_text(encoding="utf-8").splitlines():
        line = raw.strip()
        if not line:
            continue
        if line.startswith(">"):
            if name is not None:
                records.append((name, "".join(chunks).upper()))
            name = line[1:].split()[0]
            chunks = []
        else:
            chunks.append("".join(line.split()))
    if name is not None:
        records.append((name, "".join(chunks).upper()))
    return [(name, seq) for name, seq in records if seq]


def sequence_from_row(row: Mapping[str, object]) -> str:
    sequence = _string(row.get("sequence")).strip().replace(" ", "").replace("\n", "")
    if sequence:
        return sequence.upper()
    fasta_path = _string(row.get("fasta_path"))
    if fasta_path:
        records = read_fasta(fasta_path)
        if not records:
            raise ValueError(f"FASTA contains no sequence: {fasta_path}")
        return records[0][1]
    return ""


def normalize_row(raw: Mapping[str, object], manifest_path: Path) -> dict[str, str]:
    base_dir = manifest_path.parent.resolve()
    row = {str(k): _string(v) for k, v in raw.items()}
    design_id = row.get("design_id", "").strip()
    if not design_id:
        raise ValueError(f"manifest row has empty design_id: {manifest_path}")
    if not DESIGN_ID_RE.fullmatch(design_id):
        raise ValueError(
            f"invalid design_id {design_id!r}; use letters, digits, '.', '_' or '-'."
        )
    row["design_id"] = design_id
    row.setdefault("parent_design_id", "")
    row.setdefault("source_manifest", str(manifest_path.resolve()))
    for field in PATH_FIELDS:
        if field in row and row[field]:
            row[field] = _resolve(row[field], base_dir)
    if not row.get("sequence") and row.get("fasta_path"):
        row["sequence"] = sequence_from_row(row)
    row.setdefault("sequence", "")
    return row


def read_manifest(path: str | Path, *, require_sequence: bool = False) -> list[dict[str, str]]:
    manifest_path = Path(path).expanduser().resolve()
    if not manifest_path.is_file():
        raise FileNotFoundError(f"manifest not found: {manifest_path}")
    with manifest_path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames or "design_id" not in reader.fieldnames:
            raise ValueError(f"manifest must contain a design_id column: {manifest_path}")
        if len(set(reader.fieldnames)) != len(reader.fieldnames) or any(not name.strip() for name in reader.fieldnames):
            raise ValueError(f"manifest has duplicate or empty column names: {manifest_path}")
        rows = []
        for row in reader:
            if None in row or any(value is None for value in row.values()):
                raise ValueError(f"manifest line {reader.line_num}: column count differs from header")
            rows.append(normalize_row(row, manifest_path))
    if not rows:
        raise ValueError(f"manifest is empty: {manifest_path}")
    ids = [row["design_id"] for row in rows]
    duplicates = sorted({item for item in ids if ids.count(item) > 1})
    if duplicates:
        raise ValueError(f"duplicate design_id in {manifest_path}: {', '.join(duplicates)}")
    if require_sequence:
        missing = [row["design_id"] for row in rows if not row.get("sequence")]
        if missing:
            raise ValueError(f"missing sequence/fasta_path for: {', '.join(missing[:10])}")
    return rows


def _field_order(rows: Iterable[Mapping[str, object]], preferred: Iterable[str] = ()) -> list[str]:
    keys: list[str] = []
    seen: set[str] = set()
    for key in list(preferred) + list(CORE_FIELDS):
        if key not in seen:
            keys.append(key)
            seen.add(key)
    for row in rows:
        for key in row:
            if key not in seen:
                keys.append(key)
                seen.add(key)
    return keys


def write_tsv(path: str | Path, rows: Iterable[Mapping[str, object]], *, preferred: Iterable[str] = ()) -> None:
    rows_list = [dict(row) for row in rows]
    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    fields = _field_order(rows_list, preferred)
    tmp = output.with_name(f".{output.name}.tmp")
    with tmp.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows_list:
            writer.writerow({field: _string(row.get(field, "")) for field in fields})
    tmp.replace(output)


def write_csv(path: str | Path, rows: Iterable[Mapping[str, object]], *, preferred: Iterable[str] = ()) -> None:
    rows_list = [dict(row) for row in rows]
    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    fields = _field_order(rows_list, preferred)
    tmp = output.with_name(f".{output.name}.tmp")
    with tmp.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        for row in rows_list:
            writer.writerow({field: _string(row.get(field, "")) for field in fields})
    tmp.replace(output)


def write_fasta(rows: Iterable[Mapping[str, object]], path: str | Path) -> int:
    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    count = 0
    tmp = output.with_name(f".{output.name}.tmp")
    with tmp.open("w", encoding="utf-8") as handle:
        for row in rows:
            sequence = sequence_from_row(row)
            if not sequence:
                continue
            handle.write(f">{row['design_id']}\n{sequence}\n")
            count += 1
    tmp.replace(output)
    return count


def row_fingerprint(row: Mapping[str, object]) -> str:
    payload = {str(key): _string(value) for key, value in sorted(row.items())}
    return hashlib.sha256(json.dumps(payload, sort_keys=True).encode("utf-8")).hexdigest()


def metadata(row: Mapping[str, object]) -> dict[str, object]:
    raw = _string(row.get("metadata_json"))
    if not raw:
        return {}
    try:
        value = json.loads(raw)
    except json.JSONDecodeError as exc:
        raise ValueError(f"metadata_json is not valid JSON for {row.get('design_id')}") from exc
    if not isinstance(value, dict):
        raise ValueError(f"metadata_json must be an object for {row.get('design_id')}")
    return value
