#!/usr/bin/env python3
"""Dependency-light geometry filters for Level 2 structures.

These metrics are deterministic CPU fallbacks.  They do not pretend to be
Rosetta shape complementarity; a configured Rosetta batch script can add the
official Rosetta values alongside them.
"""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path
from typing import Iterable

import numpy as np


def _atoms(path: Path) -> tuple[np.ndarray, dict[str, list[int]]]:
    import gemmi
    structure = gemmi.read_structure(str(path))
    coords: list[list[float]] = []
    chains: dict[str, list[int]] = {}
    for model in structure:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    if atom.element.name == "H":
                        continue
                    index = len(coords)
                    coords.append([atom.pos.x, atom.pos.y, atom.pos.z])
                    chains.setdefault(str(chain.name), []).append(index)
        break
    return np.asarray(coords, dtype=float), chains


def _radius_gyration(coords: np.ndarray) -> float | str:
    if len(coords) == 0:
        return ""
    center = coords.mean(axis=0)
    return float(np.sqrt(np.mean(np.sum((coords - center) ** 2, axis=1))))


def _contacts(coords: np.ndarray, chains: dict[str, list[int]], cutoff: float = 5.0) -> dict[str, object]:
    names = sorted(chains)
    pairs: dict[str, object] = {}
    for left_index, left in enumerate(names):
        for right in names[left_index + 1:]:
            a = coords[chains[left]]
            b = coords[chains[right]]
            if not len(a) or not len(b):
                continue
            distance = np.sqrt(((a[:, None, :] - b[None, :, :]) ** 2).sum(axis=2))
            hits = np.argwhere(distance <= cutoff)
            if len(hits):
                left_atoms = len(set(int(item[0]) for item in hits))
                right_atoms = len(set(int(item[1]) for item in hits))
                pairs[f"{left}:{right}"] = {
                    "contact_count": int(len(hits)),
                    "left_interface_atoms": left_atoms,
                    "right_interface_atoms": right_atoms,
                    "contact_density": float(len(hits) / math.sqrt(max(1, len(a) * len(b)))),
                    "min_distance": float(distance[hits[:, 0], hits[:, 1]].min()),
                }
    return pairs


def _sasa(path: Path) -> dict[str, object]:
    """Use Biopython Shrake-Rupley when the optional dependency is present."""
    try:
        from Bio.PDB import MMCIFParser, PDBParser, ShrakeRupley
        parser = MMCIFParser(QUIET=True) if path.suffix.lower() in {".cif", ".mmcif"} else PDBParser(QUIET=True)
        structure = parser.get_structure("profold", str(path))
        ShrakeRupley().compute(structure, level="A")
        total = 0.0
        chains: dict[str, float] = {}
        for model in structure:
            for chain in model:
                value = sum(float(getattr(atom, "sasa", 0.0)) for atom in chain.get_atoms())
                chains[str(chain.id)] = chains.get(str(chain.id), 0.0) + value
                total += value
            break
        return {"sasa": total, "chain_sasa": chains}
    except Exception as exc:
        return {"sasa": "", "chain_sasa": {}, "sasa_error": f"{type(exc).__name__}: {exc}"}


def _isoelectric_point(path: Path) -> float | str:
    """Estimate protein pI from residues in the predicted structure."""
    try:
        from Bio.PDB import MMCIFParser, PDBParser
        from Bio.SeqUtils.ProtParam import ProteinAnalysis
        parser = MMCIFParser(QUIET=True) if path.suffix.lower() in {".cif", ".mmcif"} else PDBParser(QUIET=True)
        structure = parser.get_structure("profold", str(path))
        three_to_one = {
            "ALA":"A", "ARG":"R", "ASN":"N", "ASP":"D", "CYS":"C",
            "GLN":"Q", "GLU":"E", "GLY":"G", "HIS":"H", "ILE":"I",
            "LEU":"L", "LYS":"K", "MET":"M", "PHE":"F", "PRO":"P",
            "SER":"S", "THR":"T", "TRP":"W", "TYR":"Y", "VAL":"V",
        }
        sequence = []
        seen = set()
        for model in structure:
            for chain in model:
                for residue in chain:
                    key = (str(chain.id), int(residue.id[1]))
                    if key in seen or residue.id[0].strip():
                        continue
                    aa = three_to_one.get(str(residue.resname).upper())
                    if aa:
                        sequence.append(aa); seen.add(key)
            break
        return float(ProteinAnalysis("".join(sequence)).isoelectric_point()) if sequence else ""
    except Exception:
        return ""


def structure_metrics(path: str | Path, unbound_path: str | Path | None = None) -> dict[str, object]:
    structure = Path(path).expanduser().resolve()
    coords, chains = _atoms(structure)
    result: dict[str, object] = {
        "structure_path": str(structure),
        "radius_of_gyration": _radius_gyration(coords),
        "isoelectric_point": _isoelectric_point(structure),
        "chain_count": len(chains),
        "chain_atom_counts": {key: len(value) for key, value in chains.items()},
        "interface_pairs": _contacts(coords, chains),
    }
    result.update(_sasa(structure))
    if unbound_path:
        unbound = _sasa(Path(unbound_path).expanduser().resolve())
        if isinstance(result.get("sasa"), (int, float)) and isinstance(unbound.get("sasa"), (int, float)):
            result["delta_sasa"] = float(unbound["sasa"]) - float(result["sasa"])
            result["unbound_sasa"] = unbound["sasa"]
    return result


def run_geometry_filters(run_dir: str | Path) -> Path:
    run = Path(run_dir).resolve()
    manifest = run / "results" / "result_manifest.tsv"
    output = run / "metrics" / "geometry_filters.json"
    rows: list[dict[str, object]] = []
    with manifest.open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            # AF3 is the structural reference for physical filtering when a
            # Level1 result was supplied.  Fall back to each predictor's own
            # best structure for standalone Level2 runs.
            af3_path = str(row.get("af3_structure_path", "") or "").strip()
            candidates = [("af3", af3_path)] if af3_path and Path(af3_path).is_file() else []
            candidates.extend((key[:-len("_structure_path")], value) for key, value in row.items()
                              if key.endswith("_structure_path") and value and not key.startswith("af3_"))
            for model, value in candidates:
                if model in {"reference", "backbone", "structure"}:
                    continue
                path = Path(value)
                if not path.is_file():
                    continue
                try:
                    item = structure_metrics(path, row.get("unbound_structure_path") or None)
                except Exception as exc:
                    item = {"structure_path": str(path), "status": "failed", "error": f"{type(exc).__name__}: {exc}"}
                item.update({"design_id": row.get("design_id", ""), "model": model})
                rows.append(item)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(rows, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    return output
