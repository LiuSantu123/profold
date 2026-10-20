#!/usr/bin/env python3
"""Extract pLDDT and basic integrity metrics from ESMFold2 PDB outputs.

ESMFold2 deployments do not currently share a confidence JSON contract.  The
standard ESMFold convention stores pLDDT in the PDB B-factor column, so this
monitor deliberately reports only values that can be validated from the
structure itself.  It never treats a missing value as zero.
"""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path


def _pdb_metrics(path: Path) -> dict[str, object]:
    atoms = 0
    residues: set[tuple[str, str, str]] = set()
    chains: set[str] = set()
    plddt: list[float] = []
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if not (line.startswith("ATOM  ") or line.startswith("HETATM")):
            continue
        atoms += 1
        chain = line[21:22].strip() or "?"
        residue = line[17:20].strip()
        number = line[22:26].strip()
        insertion = line[26:27].strip()
        chains.add(chain)
        residues.add((chain, number, insertion))
        try:
            score = float(line[60:66].strip())
            if math.isfinite(score) and 0 <= score <= 100:
                plddt.append(score)
        except ValueError:
            pass
    if atoms == 0:
        return {"status": "failed", "error": "PDB contains no atoms"}
    if not plddt:
        return {"status": "partial", "error": "no valid B-factor pLDDT values"}
    return {
        "status": "success",
        "error": "",
        "mean_plddt": round(sum(plddt) / len(plddt), 4),
        "n_atoms": atoms,
        "n_residues": len(residues),
        "n_chains": len(chains),
    }


def scan(root: Path) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for path in sorted(root.rglob("*.pdb")):
        if not path.is_file():
            continue
        design_id = path.stem
        result = {"design_id": design_id, "structure_path": str(path.resolve())}
        try:
            result.update(_pdb_metrics(path))
        except OSError as exc:
            result.update({"status": "failed", "error": str(exc)})
        rows.append(result)
    return rows


def write_csv(rows: list[dict[str, object]], output: Path) -> None:
    """Write monitor rows using a stable column order."""
    columns = ["design_id", "status", "error", "mean_plddt", "n_atoms", "n_residues", "n_chains", "structure_path"]
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        writer.writerows({key: row.get(key, "") for key in columns} for row in rows)


def main() -> int:
    parser = argparse.ArgumentParser(description="Monitor ESMFold2 PDB outputs")
    parser.add_argument("--scan", required=True, type=Path)
    parser.add_argument("-o", "--out", required=True, type=Path)
    args = parser.parse_args()
    rows = scan(args.scan)
    if not rows:
        print(f"no PDB outputs found under {args.scan}")
        return 1
    write_csv(rows, args.out)
    ok = sum(row.get("status") == "success" for row in rows)
    print(f"completed: {ok}/{len(rows)} successful, output -> {args.out}")
    return 0 if ok == len(rows) else 2


if __name__ == "__main__":
    raise SystemExit(main())
