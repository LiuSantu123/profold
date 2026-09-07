#!/usr/bin/env python3
"""Align designed/predicted structures to their design-time backbones.

The input is deliberately a small, ordinary CSV rather than a v37-internal
manifest.  One row represents one predicted design structure and its source
backbone.  The program writes three auditable tables:

* ``alignment_summary.csv``: one row per design, including complex-level
  US-align metrics and ligand/entity metrics;
* ``chain_pairs.csv``: every backbone-chain × predicted-chain pair;
* ``motifs.csv``: exact residue/atom-pair local fits requested by the row.

US-align is used for protein/nucleic-acid polymer alignments.  Non-polymer
ligands are not silently passed to US-align: their heavy atoms are matched by
residue/atom name and scored after the polymer-derived superposition.  This
keeps the complex-level US-align score and ligand RMSD scientifically
separate while still supporting protein–DNA/RNA–small-molecule complexes.

No prediction, MSA, or source-file mutation happens here.
"""

from __future__ import annotations

import argparse
import copy
import csv
import difflib
import json
import math
import os
import re
import subprocess
import sys
import tempfile
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

try:
    import gemmi
except ImportError as exc:  # pragma: no cover - exercised only on a bad host
    raise SystemExit("structure_align requires gemmi; use the v37/base Python environment") from exc

import numpy as np


DEFAULT_USALIGN = "/xcfhome/yhliu/002_software/009_TMalign/USalign"
WATER_NAMES = {"HOH", "WAT", "DOD", "H2O"}
NUCLEIC_NAMES = {
    "A", "C", "G", "U", "I", "N", "T",
    "DA", "DC", "DG", "DT", "DU", "DI",
    "RA", "RC", "RG", "RU", "RI",
    "PSU", "5MC", "OMC", "OMG", "1MA", "2MG", "M2G", "7MG",
}
AA3_TO_1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    "MSE": "M", "SEC": "U", "PYL": "O",
}
NUC_TO_1 = {
    "A": "A", "DA": "A", "RA": "A", "C": "C", "DC": "C", "RC": "C",
    "G": "G", "DG": "G", "RG": "G", "U": "U", "DU": "U", "RU": "U",
    "T": "T", "DT": "T", "I": "I", "DI": "I", "RI": "I",
}


def _safe_float(value: object) -> float | None:
    try:
        number = float(value)  # type: ignore[arg-type]
    except (TypeError, ValueError):
        return None
    return number if math.isfinite(number) else None


def _fmt_number(value: float | None, digits: int = 5) -> object:
    if value is None or not math.isfinite(value):
        return ""
    return round(float(value), digits)


def _coord(atom: Any) -> np.ndarray:
    return np.asarray([atom.pos.x, atom.pos.y, atom.pos.z], dtype=float)


def _is_hydrogen(atom: Any) -> bool:
    try:
        return bool(atom.is_hydrogen())
    except AttributeError:
        return str(getattr(atom, "element", "")).upper() == "H"


def _residue_kind(residue: Any) -> str:
    if residue.entity_type != gemmi.EntityType.Polymer:
        return "ligand"
    name = residue.name.strip().upper()
    if name in NUCLEIC_NAMES:
        return "nucleic_acid"
    return "protein"


def _atom_element(atom: Any) -> str:
    element = getattr(atom, "element", "")
    try:
        return str(element.name).upper()
    except AttributeError:
        return str(element).upper()


def _preferred_atom(residue: Any, kind: str) -> Any | None:
    atoms = {str(atom.name).strip().upper(): atom for atom in residue if not _is_hydrogen(atom)}
    if kind == "protein":
        names = ("CA",)
    elif kind == "nucleic_acid":
        names = ("C4'", "C4*", "C1'", "C1*", "C1")
    else:
        names = ()
    for name in names:
        if name in atoms:
            return atoms[name]
    for atom in residue:
        if not _is_hydrogen(atom):
            return atom
    return None


@dataclass
class AtomInfo:
    chain_id: str
    residue_name: str
    residue_num: int
    icode: str
    atom_name: str
    element: str
    coord: np.ndarray

    @property
    def residue_key(self) -> str:
        return f"{self.residue_num}{self.icode}"


@dataclass
class ResidueInfo:
    chain_id: str
    name: str
    residue_num: int
    icode: str
    kind: str
    token_coord: np.ndarray | None
    atom_coords: dict[str, np.ndarray] = field(default_factory=dict)

    @property
    def key(self) -> str:
        return f"{self.residue_num}{self.icode}"

    @property
    def one_letter(self) -> str:
        if self.kind == "nucleic_acid":
            return NUC_TO_1.get(self.name, "N")
        return AA3_TO_1.get(self.name, "X")


@dataclass
class ChainInfo:
    chain_id: str
    residues: list[ResidueInfo]
    ligand_atoms: list[AtomInfo]
    kind: str

    @property
    def polymer_residues(self) -> list[ResidueInfo]:
        return [residue for residue in self.residues if residue.kind in {"protein", "nucleic_acid"}]

    @property
    def polymer_kind(self) -> str:
        kinds = {residue.kind for residue in self.residues if residue.kind in {"protein", "nucleic_acid"}}
        if len(kinds) == 1:
            return next(iter(kinds))
        if len(kinds) > 1:
            return "mixed_polymer"
        return ""

    @property
    def sequence(self) -> str:
        return "".join(residue.one_letter for residue in self.polymer_residues)


@dataclass
class StructureInfo:
    path: Path
    chains: list[ChainInfo]

    @property
    def chain_by_id(self) -> dict[str, ChainInfo]:
        return {chain.chain_id: chain for chain in self.chains}

    @property
    def polymer_chains(self) -> list[ChainInfo]:
        return [chain for chain in self.chains if chain.polymer_residues]

    @property
    def ligand_chains(self) -> list[ChainInfo]:
        return [chain for chain in self.chains if chain.ligand_atoms]


@dataclass
class USAlignResult:
    status: str
    error: str = ""
    aligned_length: int | None = None
    rmsd: float | None = None
    seq_identity: float | None = None
    tm_score_structure1: float | None = None
    tm_score_structure2: float | None = None
    raw: str = ""
    rotation: np.ndarray | None = None
    translation: np.ndarray | None = None


def _parse_structure(path: Path) -> StructureInfo:
    if not path.is_file():
        raise FileNotFoundError(path)
    structure = gemmi.read_structure(str(path))
    if len(structure) == 0:
        raise ValueError(f"structure has no model: {path}")
    chains: list[ChainInfo] = []
    for chain in structure[0]:
        chain_id = (chain.name or "_").strip() or "_"
        residues: list[ResidueInfo] = []
        ligand_atoms: list[AtomInfo] = []
        for residue in chain:
            name = residue.name.strip().upper()
            if name in WATER_NAMES:
                continue
            kind = _residue_kind(residue)
            seq_num = int(residue.seqid.num)
            icode = str(residue.seqid.icode or "").strip()
            atoms: dict[str, np.ndarray] = {}
            for atom in residue:
                if _is_hydrogen(atom):
                    continue
                atoms[str(atom.name).strip().upper()] = _coord(atom)
                if kind == "ligand":
                    ligand_atoms.append(
                        AtomInfo(
                            chain_id=chain_id,
                            residue_name=name,
                            residue_num=seq_num,
                            icode=icode,
                            atom_name=str(atom.name).strip().upper(),
                            element=_atom_element(atom),
                            coord=_coord(atom),
                        )
                    )
            if kind == "ligand" and not atoms:
                continue
            preferred = _preferred_atom(residue, kind)
            residues.append(
                ResidueInfo(
                    chain_id=chain_id,
                    name=name,
                    residue_num=seq_num,
                    icode=icode,
                    kind=kind,
                    token_coord=_coord(preferred) if preferred is not None else None,
                    atom_coords=atoms,
                )
            )
        if not residues and not ligand_atoms:
            continue
        polymer_kinds = {residue.kind for residue in residues if residue.kind in {"protein", "nucleic_acid"}}
        if not polymer_kinds:
            chain_kind = "ligand"
        elif len(polymer_kinds) == 1 and not ligand_atoms:
            chain_kind = next(iter(polymer_kinds))
        else:
            chain_kind = "mixed"
        chains.append(ChainInfo(chain_id, residues, ligand_atoms, chain_kind))
    if not chains:
        raise ValueError(f"structure has no usable chains/entities: {path}")
    return StructureInfo(path, chains)


def _sequence_similarity(first: ChainInfo, second: ChainInfo) -> float:
    if not first.sequence or not second.sequence:
        return 0.0
    ratio = difflib.SequenceMatcher(None, first.sequence, second.sequence, autojunk=False).ratio()
    length_ratio = min(len(first.sequence), len(second.sequence)) / max(len(first.sequence), len(second.sequence))
    kind_bonus = 0.12 if first.polymer_kind == second.polymer_kind else -0.25
    return max(-1.0, 0.72 * ratio + 0.28 * length_ratio + kind_bonus)


def _ligand_similarity(first: ChainInfo, second: ChainInfo) -> float:
    if not first.ligand_atoms or not second.ligand_atoms:
        return -1.0
    first_names = {atom.residue_name for atom in first.ligand_atoms}
    second_names = {atom.residue_name for atom in second.ligand_atoms}
    name_score = 1.0 if first_names & second_names else 0.0
    length_score = min(len(first.ligand_atoms), len(second.ligand_atoms)) / max(len(first.ligand_atoms), len(second.ligand_atoms))
    return 0.75 * name_score + 0.25 * length_score


def _chain_pair_score(backbone: ChainInfo, predicted: ChainInfo) -> float:
    if backbone.polymer_residues and predicted.polymer_residues:
        return _sequence_similarity(backbone, predicted)
    if backbone.ligand_atoms and predicted.ligand_atoms:
        return _ligand_similarity(backbone, predicted)
    return -1.0


def auto_chain_mapping(backbone: StructureInfo, predicted: StructureInfo) -> dict[str, str]:
    """Greedily choose a transparent one-to-one chain map.

    All chain pairs are still reported separately.  This map only marks the
    most plausible pair for motif defaults and ligand superposition.
    """

    candidates: list[tuple[float, str, str]] = []
    for ref in backbone.chains:
        for pred in predicted.chains:
            score = _chain_pair_score(ref, pred)
            if score >= 0:
                candidates.append((score, ref.chain_id, pred.chain_id))
    mapping: dict[str, str] = {}
    used_pred: set[str] = set()
    for _score, ref_id, pred_id in sorted(candidates, key=lambda item: (-item[0], item[1], item[2])):
        if ref_id in mapping or pred_id in used_pred:
            continue
        mapping[ref_id] = pred_id
        used_pred.add(pred_id)
    return mapping


def _parse_chain_map(value: object, base_dir: Path) -> dict[str, str]:
    if value in (None, ""):
        return {}
    text = str(value).strip()
    candidate = Path(text).expanduser()
    if text.startswith("@"):
        candidate = Path(text[1:]).expanduser()
    if not candidate.is_absolute():
        candidate = (base_dir / candidate).resolve()
    if text.startswith("@") or candidate.is_file():
        text = candidate.read_text(encoding="utf-8")
    try:
        payload = json.loads(text)
        if isinstance(payload, dict):
            return {str(key).strip(): str(val).strip() for key, val in payload.items() if str(key).strip() and str(val).strip()}
        if isinstance(payload, list):
            output: dict[str, str] = {}
            for item in payload:
                if isinstance(item, dict) and item.get("backbone_chain") and item.get("predicted_chain"):
                    output[str(item["backbone_chain"])] = str(item["predicted_chain"])
            if output:
                return output
    except json.JSONDecodeError:
        pass
    output = {}
    for item in re.split(r"[,;\s]+", text):
        item = item.strip()
        if not item:
            continue
        match = re.match(r"^([^:=\- >]+)\s*(?:->|:|=)\s*([^:=\- >]+)$", item)
        if not match:
            raise ValueError(f"invalid chain_map item {item!r}; use A:B,C:D or JSON")
        output[match.group(1)] = match.group(2)
    return output


def _explicit_or_auto_map(row: Mapping[str, str], backbone: StructureInfo, predicted: StructureInfo, base_dir: Path) -> tuple[dict[str, str], str]:
    value = row.get("backbone_to_predicted_chain_map") or row.get("chain_map") or row.get("backbone_chain_map")
    if value:
        mapping = _parse_chain_map(value, base_dir)
        for ref_id, pred_id in mapping.items():
            if ref_id not in backbone.chain_by_id:
                raise ValueError(f"chain_map references missing backbone chain {ref_id}")
            if pred_id not in predicted.chain_by_id:
                raise ValueError(f"chain_map references missing predicted chain {pred_id}")
        if len(set(mapping.values())) != len(mapping):
            raise ValueError("chain_map maps more than one backbone chain to a predicted chain")
        return mapping, "explicit"
    return auto_chain_mapping(backbone, predicted), "auto"


def _parse_usalign(text: str, returncode: int, error: str = "") -> USAlignResult:
    aligned = re.search(r"Aligned length=\s*(\d+),\s*RMSD=\s*([0-9.eE+-]+),\s*Seq_ID=.*?=\s*([0-9.eE+-]+)", text)
    tm_scores = [float(match) for match in re.findall(r"TM-score=\s*([0-9.eE+-]+)", text)]
    if returncode != 0:
        message = error or f"USalign returncode={returncode}"
        return USAlignResult("failed", message, raw=text)
    if not aligned:
        return USAlignResult("failed", "USalign output has no aligned-length line", raw=text)
    return USAlignResult(
        status="success",
        aligned_length=int(aligned.group(1)),
        rmsd=float(aligned.group(2)),
        seq_identity=float(aligned.group(3)),
        tm_score_structure1=tm_scores[0] if tm_scores else None,
        tm_score_structure2=tm_scores[1] if len(tm_scores) > 1 else (tm_scores[0] if tm_scores else None),
        raw=text,
    )


def _parse_rotation_matrix(path: Path) -> tuple[np.ndarray | None, np.ndarray | None]:
    if not path.is_file():
        return None, None
    rows: list[list[float]] = []
    try:
        for raw in path.read_text(encoding="utf-8", errors="replace").splitlines():
            fields = raw.split()
            if len(fields) != 5 or not fields[0].isdigit():
                continue
            try:
                rows.append([float(item) for item in fields[1:]])
            except ValueError:
                continue
    except OSError:
        return None, None
    if len(rows) < 3:
        return None, None
    values = np.asarray(rows[:3], dtype=float)
    return values[:, 1:4], values[:, 0]


def run_usalign(
    predicted_path: Path,
    backbone_path: Path,
    usalign_path: Path,
    *,
    multi_chain: bool,
    molecule: str = "auto",
    timeout: int = 300,
    matrix_path: Path | None = None,
) -> USAlignResult:
    if not usalign_path.is_file() or not os.access(usalign_path, os.X_OK):
        return USAlignResult("failed", f"USalign is not executable: {usalign_path}")
    command = [str(usalign_path), str(predicted_path), str(backbone_path)]
    if molecule == "nucleic_acid":
        command.extend(["-mol", "RNA"])
    if multi_chain:
        command.extend(["-mm", "1", "-ter", "1"])
    else:
        command.extend(["-ter", "1"])
    if matrix_path is not None:
        command.extend(["-m", str(matrix_path)])
    try:
        completed = subprocess.run(
            command,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            check=False,
            timeout=timeout,
        )
    except subprocess.TimeoutExpired as exc:
        raw = (exc.stdout or "") if isinstance(exc.stdout, str) else ""
        return USAlignResult("failed", f"USalign timeout after {timeout}s", raw=raw)
    except OSError as exc:
        return USAlignResult("failed", repr(exc))
    result = _parse_usalign(completed.stdout or "", completed.returncode)
    if matrix_path is not None:
        result.rotation, result.translation = _parse_rotation_matrix(matrix_path)
    return result


def _write_chain_pdb(source: Path, chain_id: str, output: Path) -> None:
    structure = gemmi.read_structure(str(source))
    if len(structure) == 0:
        raise ValueError(f"no model in {source}")
    chain = structure[0].find_chain(chain_id)
    if chain is None:
        raise ValueError(f"chain {chain_id!r} not found in {source}")
    selected = gemmi.Structure()
    selected.cell = structure.cell
    selected.add_model(gemmi.Model("1"))
    selected[0].add_chain(copy.deepcopy(chain))
    output.parent.mkdir(parents=True, exist_ok=True)
    selected.write_pdb(str(output))


def _match_ligand_atoms(first: ChainInfo, second: ChainInfo) -> tuple[np.ndarray, np.ndarray]:
    """Return predicted/reference atom coordinates with conservative matching."""

    def grouped(chain: ChainInfo) -> dict[tuple[str, str, str], list[AtomInfo]]:
        result: dict[tuple[str, str, str], list[AtomInfo]] = {}
        for atom in chain.ligand_atoms:
            key = (atom.residue_name, atom.atom_name, atom.element)
            result.setdefault(key, []).append(atom)
        return result

    first_groups = grouped(first)
    second_groups = grouped(second)
    first_coords: list[np.ndarray] = []
    second_coords: list[np.ndarray] = []
    for key in sorted(set(first_groups) & set(second_groups)):
        left = sorted(first_groups[key], key=lambda item: (item.residue_num, item.icode))
        right = sorted(second_groups[key], key=lambda item: (item.residue_num, item.icode))
        for a, b in zip(left, right):
            first_coords.append(a.coord)
            second_coords.append(b.coord)
    if first_coords:
        return np.vstack(first_coords), np.vstack(second_coords)
    # Some prediction writers change atom names but preserve residue and
    # element.  A second, still conservative, fallback handles that case.
    by_residue: dict[tuple[str, str], list[AtomInfo]] = {}
    by_residue_ref: dict[tuple[str, str], list[AtomInfo]] = {}
    for atom in first.ligand_atoms:
        by_residue.setdefault((atom.residue_name, atom.element), []).append(atom)
    for atom in second.ligand_atoms:
        by_residue_ref.setdefault((atom.residue_name, atom.element), []).append(atom)
    for key in sorted(set(by_residue) & set(by_residue_ref)):
        left = sorted(by_residue[key], key=lambda item: (item.residue_num, item.icode, item.atom_name))
        right = sorted(by_residue_ref[key], key=lambda item: (item.residue_num, item.icode, item.atom_name))
        for a, b in zip(left, right):
            first_coords.append(a.coord)
            second_coords.append(b.coord)
    if not first_coords:
        return np.empty((0, 3)), np.empty((0, 3))
    return np.vstack(first_coords), np.vstack(second_coords)


def _kabsch(predicted: np.ndarray, reference: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    if predicted.shape != reference.shape or predicted.ndim != 2 or predicted.shape[1] != 3:
        raise ValueError("Kabsch coordinate arrays must have identical shape (N, 3)")
    if len(predicted) == 0:
        return np.empty((0, 3)), np.eye(3), np.zeros(3)
    pred_center = predicted.mean(axis=0)
    ref_center = reference.mean(axis=0)
    if len(predicted) < 2:
        rotation = np.eye(3)
    else:
        covariance = (predicted - pred_center).T @ (reference - ref_center)
        left, _singular, right_t = np.linalg.svd(covariance)
        rotation = right_t.T @ left.T
        if np.linalg.det(rotation) < 0:
            right_t[-1, :] *= -1
            rotation = right_t.T @ left.T
    transformed = (predicted - pred_center) @ rotation + ref_center
    return transformed, rotation, ref_center - pred_center @ rotation


def _tm_score_from_distances(distances: np.ndarray, length: int) -> float | None:
    if length <= 0:
        return None
    if length <= 21:
        d0 = 0.5
    else:
        d0 = max(0.5, 1.24 * ((length - 15) ** (1.0 / 3.0)) - 1.8)
    return float(np.mean(1.0 / (1.0 + (distances / d0) ** 2)))


def _coordinate_metrics(predicted: np.ndarray, reference: np.ndarray) -> dict[str, object]:
    if len(predicted) == 0:
        return {"status": "missing", "error": "no matched coordinates", "n_pairs": 0}
    transformed, rotation, translation = _kabsch(predicted, reference)
    distances = np.linalg.norm(transformed - reference, axis=1)
    return {
        "status": "success",
        "error": "",
        "n_pairs": int(len(distances)),
        "rmsd": _fmt_number(float(np.sqrt(np.mean(distances**2)))),
        "tm_score": _fmt_number(_tm_score_from_distances(distances, len(distances))),
        "raw_rmsd": _fmt_number(float(np.sqrt(np.mean((predicted - reference) ** 2)))),
        "rotation": rotation,
        "translation": translation,
    }


def _apply_usalign_transform(coords: np.ndarray, result: USAlignResult) -> np.ndarray:
    if result.rotation is None or result.translation is None:
        raise ValueError("USalign did not provide a rotation matrix")
    return coords @ result.rotation.T + result.translation


def _residue_lookup(chain: ChainInfo) -> dict[str, ResidueInfo]:
    return {residue.key: residue for residue in chain.polymer_residues}


def _parse_selector(value: object) -> list[str]:
    if value in (None, ""):
        return []
    if isinstance(value, list):
        values = value
    else:
        text = str(value).strip()
        try:
            parsed = json.loads(text)
            values = parsed if isinstance(parsed, list) else [parsed]
        except json.JSONDecodeError:
            values = re.split(r"[,;\s]+", text)
    output: list[str] = []
    for item in values:
        token = str(item).strip()
        if not token:
            continue
        match = re.fullmatch(r"(-?\d+)([A-Za-z]?)\s*-\s*(-?\d+)([A-Za-z]?)", token)
        if match:
            start, end = int(match.group(1)), int(match.group(3))
            suffix_start, suffix_end = match.group(2), match.group(4)
            if suffix_start or suffix_end:
                if start == end:
                    output.append(f"{start}{suffix_start or suffix_end}")
                else:
                    output.extend(str(number) for number in range(start, end + (1 if end >= start else -1), 1 if end >= start else -1))
            else:
                step = 1 if end >= start else -1
                output.extend(str(number) for number in range(start, end + step, step))
        else:
            output.append(token)
    return output


def _load_motif_specs(row: Mapping[str, str], base_dir: Path) -> list[dict[str, Any]]:
    raw = row.get("motif_json") or row.get("motif_spec") or row.get("motif")
    if raw:
        text = str(raw).strip()
        candidate = Path(text[1:] if text.startswith("@") else text).expanduser()
        if not candidate.is_absolute():
            candidate = (base_dir / candidate).resolve()
        if text.startswith("@") or candidate.is_file():
            text = candidate.read_text(encoding="utf-8")
        try:
            payload = json.loads(text)
        except json.JSONDecodeError:
            payload = None
        if isinstance(payload, dict):
            payload = payload.get("motifs", [payload])
        if isinstance(payload, list):
            return [dict(item) for item in payload if isinstance(item, dict)]
        # Compact syntax: name:chain:residues;name2:chain:residues
        specs: list[dict[str, Any]] = []
        for item in text.split(";"):
            fields = [part.strip() for part in item.split(":")]
            if len(fields) >= 3:
                specs.append({"name": fields[0], "chain": fields[1], "residues": fields[2]})
        if specs:
            return specs
    if row.get("motif_name") or row.get("motif_residues") or row.get("motif_chain"):
        return [{
            "name": row.get("motif_name") or "motif",
            "chain": row.get("motif_chain", ""),
            "residues": row.get("motif_residues", ""),
            "atom": row.get("motif_atom", ""),
        }]
    return []


def _residues_from_spec(chain: ChainInfo, spec: Mapping[str, Any], side: str) -> list[ResidueInfo]:
    lookup = _residue_lookup(chain)
    key = f"{side}_residues"
    selectors = spec.get(key)
    if selectors in (None, ""):
        selectors = spec.get("residues")
    if selectors not in (None, ""):
        result = []
        for selector in _parse_selector(selectors):
            residue = lookup.get(str(selector))
            if residue is not None:
                result.append(residue)
        return result
    positions = spec.get(f"{side}_positions") or spec.get("positions")
    if positions in (None, ""):
        return []
    result = []
    polymer = chain.polymer_residues
    for token in _parse_selector(positions):
        try:
            index = int(token) - 1
        except ValueError:
            continue
        if 0 <= index < len(polymer):
            result.append(polymer[index])
    return result


def _motif_coords(
    backbone: ChainInfo,
    predicted: ChainInfo,
    spec: Mapping[str, Any],
) -> tuple[np.ndarray, np.ndarray, int]:
    reference_residues = _residues_from_spec(backbone, spec, "backbone")
    predicted_residues = _residues_from_spec(predicted, spec, "predicted")
    if not reference_residues or not predicted_residues:
        return np.empty((0, 3)), np.empty((0, 3)), 0
    # If both sides used the same selector, order is already meaningful.  If
    # only one side was given, pair in order and make the truncation explicit.
    n = min(len(reference_residues), len(predicted_residues))
    reference_residues = reference_residues[:n]
    predicted_residues = predicted_residues[:n]
    atom_selector = str(spec.get("atom") or spec.get("atom_name") or "auto").strip().upper()
    predicted_coords: list[np.ndarray] = []
    reference_coords: list[np.ndarray] = []
    for ref_residue, pred_residue in zip(reference_residues, predicted_residues):
        if atom_selector in {"ALL", "HEAVY", "ALL_HEAVY"}:
            common = sorted(set(ref_residue.atom_coords) & set(pred_residue.atom_coords))
            for atom_name in common:
                predicted_coords.append(pred_residue.atom_coords[atom_name])
                reference_coords.append(ref_residue.atom_coords[atom_name])
            continue
        if atom_selector == "AUTO":
            atom_selector_for_pair = "CA" if ref_residue.kind == "protein" else "C4'"
            candidates = [atom_selector_for_pair, "C4*", "C1'", "C1*", "C1"]
        else:
            candidates = [atom_selector]
        for atom_name in candidates:
            if atom_name in ref_residue.atom_coords and atom_name in pred_residue.atom_coords:
                predicted_coords.append(pred_residue.atom_coords[atom_name])
                reference_coords.append(ref_residue.atom_coords[atom_name])
                break
    if not predicted_coords:
        return np.empty((0, 3)), np.empty((0, 3)), n
    return np.vstack(predicted_coords), np.vstack(reference_coords), n


def _summary_base(row: Mapping[str, str], backbone: StructureInfo, predicted: StructureInfo) -> dict[str, object]:
    return {
        "design_id": row.get("design_id", ""),
        "sequence_id": row.get("sequence_id") or row.get("design_sequence_id") or row.get("design_id", ""),
        "backbone_structure_path": str(backbone.path),
        "predicted_structure_path": str(predicted.path),
        "status": "success",
        "error": "",
        "backbone_chain_count": len(backbone.chains),
        "predicted_chain_count": len(predicted.chains),
        "backbone_polymer_chain_count": len(backbone.polymer_chains),
        "predicted_polymer_chain_count": len(predicted.polymer_chains),
    }


def align_one_row(
    row: Mapping[str, str],
    *,
    input_base: Path,
    output_dir: Path,
    usalign_path: Path,
    timeout: int = 300,
    keep_raw: bool = False,
) -> tuple[dict[str, object], list[dict[str, object]], list[dict[str, object]]]:
    """Align one CSV row and return summary, chain-pair, and motif rows."""

    design_id = row.get("design_id", "")
    backbone_value = row.get("backbone_structure_path") or row.get("reference_structure_path") or row.get("backbone_path")
    predicted_value = row.get("predicted_structure_path") or row.get("structure_path") or row.get("prediction_structure_path")
    if not backbone_value or not predicted_value:
        summary = {"design_id": design_id, "status": "failed", "error": "missing backbone_structure_path or predicted_structure_path"}
        return summary, [], []
    backbone_path = Path(backbone_value).expanduser()
    predicted_path = Path(predicted_value).expanduser()
    if not backbone_path.is_absolute():
        backbone_path = (input_base / backbone_path).resolve()
    if not predicted_path.is_absolute():
        predicted_path = (input_base / predicted_path).resolve()
    raw_dir = output_dir / "raw" / design_id
    raw_dir.mkdir(parents=True, exist_ok=True)
    try:
        backbone = _parse_structure(backbone_path)
        predicted = _parse_structure(predicted_path)
        summary = _summary_base({**dict(row), "design_id": design_id}, backbone, predicted)
        mapping, mapping_method = _explicit_or_auto_map(row, backbone, predicted, input_base)
        summary["chain_mapping_method"] = mapping_method
        summary["chain_mapping_json"] = json.dumps(mapping, sort_keys=True)
        summary["backbone_chain_ids"] = ",".join(chain.chain_id for chain in backbone.chains)
        summary["predicted_chain_ids"] = ",".join(chain.chain_id for chain in predicted.chains)

        with tempfile.TemporaryDirectory(prefix=f"v37_align_{design_id}_", dir=str(output_dir)) as tmp_name:
            temp_dir = Path(tmp_name)
            matrix_path = temp_dir / "overall.matrix"
            overall: USAlignResult | None = None
            if backbone.polymer_chains and predicted.polymer_chains:
                overall = run_usalign(
                    predicted_path,
                    backbone_path,
                    usalign_path,
                    multi_chain=len(backbone.polymer_chains) > 1 or len(predicted.polymer_chains) > 1,
                    timeout=timeout,
                    matrix_path=matrix_path,
                )
                if keep_raw:
                    (raw_dir / "overall.usalign.txt").write_text(overall.raw, encoding="utf-8")
            if overall is None:
                summary.update({
                    "overall_usalign_status": "missing",
                    "overall_usalign_error": "no polymer chains in both structures",
                })
            else:
                summary.update({
                    "overall_usalign_status": overall.status,
                    "overall_usalign_error": overall.error,
                    "overall_aligned_length": overall.aligned_length,
                    "overall_usalign_rmsd": _fmt_number(overall.rmsd),
                    "overall_seq_identity": _fmt_number(overall.seq_identity),
                    "overall_tmscore_norm_predicted": _fmt_number(overall.tm_score_structure1),
                    "overall_tmscore_norm_backbone": _fmt_number(overall.tm_score_structure2),
                })

            chain_rows: list[dict[str, object]] = []
            chain_files: dict[tuple[str, str], tuple[Path, Path]] = {}
            for ref_chain in backbone.chains:
                for pred_chain in predicted.chains:
                    pair: dict[str, object] = {
                        "design_id": design_id,
                        "backbone_chain": ref_chain.chain_id,
                        "predicted_chain": pred_chain.chain_id,
                        "backbone_entity_kind": ref_chain.kind,
                        "predicted_entity_kind": pred_chain.kind,
                        "backbone_polymer_kind": ref_chain.polymer_kind,
                        "predicted_polymer_kind": pred_chain.polymer_kind,
                        "is_selected_mapping": int(mapping.get(ref_chain.chain_id) == pred_chain.chain_id),
                        "status": "missing",
                        "error": "",
                        "alignment_method": "",
                    }
                    if ref_chain.polymer_residues and pred_chain.polymer_residues:
                        if ref_chain.polymer_kind != pred_chain.polymer_kind and "mixed_polymer" not in {ref_chain.polymer_kind, pred_chain.polymer_kind}:
                            pair.update({"status": "incompatible", "error": "protein/nucleic-acid chain kinds differ"})
                        else:
                            key = (ref_chain.chain_id, pred_chain.chain_id)
                            if key not in chain_files:
                                ref_file = temp_dir / f"backbone_{len(chain_files):04d}.pdb"
                                pred_file = temp_dir / f"predicted_{len(chain_files):04d}.pdb"
                                _write_chain_pdb(backbone_path, ref_chain.chain_id, ref_file)
                                _write_chain_pdb(predicted_path, pred_chain.chain_id, pred_file)
                                chain_files[key] = (ref_file, pred_file)
                            ref_file, pred_file = chain_files[key]
                            molecule = "nucleic_acid" if ref_chain.polymer_kind == "nucleic_acid" else "auto"
                            result = run_usalign(pred_file, ref_file, usalign_path, multi_chain=False, molecule=molecule, timeout=timeout)
                            pair.update({
                                "status": result.status,
                                "error": result.error,
                                "alignment_method": "USalign",
                                "aligned_length": result.aligned_length,
                                "rmsd": _fmt_number(result.rmsd),
                                "seq_identity": _fmt_number(result.seq_identity),
                                "tm_score_norm_predicted": _fmt_number(result.tm_score_structure1),
                                "tm_score_norm_backbone": _fmt_number(result.tm_score_structure2),
                            })
                            if keep_raw:
                                pair_raw = raw_dir / f"chain_{ref_chain.chain_id}_to_{pred_chain.chain_id}.usalign.txt"
                                pair_raw.write_text(result.raw, encoding="utf-8")
                    elif ref_chain.ligand_atoms and pred_chain.ligand_atoms:
                        pred_coords, ref_coords = _match_ligand_atoms(pred_chain, ref_chain)
                        metrics = _coordinate_metrics(pred_coords, ref_coords)
                        pair.update({
                            "status": metrics["status"],
                            "error": metrics["error"],
                            "alignment_method": "ligand_atom_kabsch",
                            "n_pairs": metrics["n_pairs"],
                            "rmsd": metrics.get("rmsd", ""),
                            "tm_score_norm_backbone": metrics.get("tm_score", ""),
                            "tm_score_method": "coordinate_fit",
                        })
                    else:
                        pair.update({"status": "incompatible", "error": "polymer/non-polymer entity types differ"})
                    chain_rows.append(pair)

            ligand_pred: list[np.ndarray] = []
            ligand_ref: list[np.ndarray] = []
            for ref_id, pred_id in mapping.items():
                ref_chain = backbone.chain_by_id[ref_id]
                pred_chain = predicted.chain_by_id[pred_id]
                if not ref_chain.ligand_atoms or not pred_chain.ligand_atoms:
                    continue
                pred_coords, ref_coords = _match_ligand_atoms(pred_chain, ref_chain)
                if len(pred_coords) == 0:
                    continue
                if overall is not None and overall.status == "success" and overall.rotation is not None:
                    ligand_pred.append(_apply_usalign_transform(pred_coords, overall))
                    ligand_ref.append(ref_coords)
                else:
                    ligand_pred.append(pred_coords)
                    ligand_ref.append(ref_coords)
            if ligand_pred:
                pred_atoms = np.vstack(ligand_pred)
                ref_atoms = np.vstack(ligand_ref)
                if overall is not None and overall.status == "success" and overall.rotation is not None:
                    distances = np.linalg.norm(pred_atoms - ref_atoms, axis=1)
                    summary.update({
                        "overall_ligand_status": "success",
                        "overall_ligand_rmsd": _fmt_number(float(np.sqrt(np.mean(distances**2)))),
                        "overall_ligand_atoms": int(len(distances)),
                        "overall_ligand_fit_method": "USalign_polymer_transform",
                    })
                else:
                    metrics = _coordinate_metrics(pred_atoms, ref_atoms)
                    summary.update({
                        "overall_ligand_status": metrics["status"],
                        "overall_ligand_rmsd": metrics.get("rmsd", ""),
                        "overall_ligand_atoms": metrics["n_pairs"],
                        "overall_ligand_fit_method": "ligand_atom_kabsch",
                    })
            else:
                summary.update({"overall_ligand_status": "missing", "overall_ligand_rmsd": "", "overall_ligand_atoms": 0})

            motif_rows: list[dict[str, object]] = []
            for index, spec in enumerate(_load_motif_specs(row, input_base), start=1):
                name = str(spec.get("name") or f"motif_{index}")
                ref_id = str(spec.get("backbone_chain") or spec.get("chain") or "")
                pred_id = str(spec.get("predicted_chain") or "")
                if not pred_id and ref_id:
                    pred_id = mapping.get(ref_id, "")
                motif: dict[str, object] = {
                    "design_id": design_id,
                    "motif_name": name,
                    "backbone_chain": ref_id,
                    "predicted_chain": pred_id,
                    "status": "failed",
                    "error": "",
                }
                if ref_id not in backbone.chain_by_id or pred_id not in predicted.chain_by_id:
                    motif["error"] = "motif chain not found or could not be mapped"
                    motif_rows.append(motif)
                    continue
                pred_coords, ref_coords, selected_residues = _motif_coords(backbone.chain_by_id[ref_id], predicted.chain_by_id[pred_id], spec)
                metrics = _coordinate_metrics(pred_coords, ref_coords)
                motif.update({
                    "status": metrics["status"],
                    "error": metrics["error"],
                    "n_residues_requested": selected_residues,
                    "n_pairs": metrics["n_pairs"],
                    "rmsd": metrics.get("rmsd", ""),
                    "tm_score": metrics.get("tm_score", ""),
                    "raw_rmsd": metrics.get("raw_rmsd", ""),
                    "method": "paired_atom_kabsch",
                })
                motif_rows.append(motif)
            summary["motif_count"] = len(motif_rows)
            summary["motif_success_count"] = sum(item.get("status") == "success" for item in motif_rows)
            summary["chain_pair_count"] = len(chain_rows)
            summary["chain_pair_success_count"] = sum(item.get("status") == "success" for item in chain_rows)
            if summary.get("overall_usalign_status") not in {"success", "missing"} and not summary.get("overall_ligand_atoms"):
                summary["status"] = "partial"
            elif summary.get("overall_usalign_status") == "missing" and not summary.get("overall_ligand_atoms"):
                summary["status"] = "partial"
            return summary, chain_rows, motif_rows
    except Exception as exc:
        return {
            "design_id": design_id,
            "backbone_structure_path": str(backbone_path),
            "predicted_structure_path": str(predicted_path),
            "status": "failed",
            "error": f"{type(exc).__name__}: {exc}",
        }, [], []


def _read_input_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        if not reader.fieldnames or "design_id" not in reader.fieldnames:
            raise ValueError("input CSV must contain design_id")
        rows = [{str(key): (value or "") for key, value in row.items()} for row in reader]
    if not rows:
        raise ValueError(f"input CSV is empty: {path}")
    ids = [row["design_id"] for row in rows]
    duplicates = sorted({item for item in ids if ids.count(item) > 1})
    if duplicates:
        raise ValueError(f"duplicate design_id: {', '.join(duplicates)}")
    return rows


def _write_rows(path: Path, rows: Iterable[Mapping[str, object]], preferred: Sequence[str] = ()) -> None:
    values = [dict(row) for row in rows]
    fields: list[str] = []
    seen: set[str] = set()
    for key in list(preferred) + [key for row in values for key in row]:
        if key not in seen:
            fields.append(key)
            seen.add(key)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(f".{path.name}.tmp")
    with tmp.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        for row in values:
            writer.writerow({field: row.get(field, "") for field in fields})
    tmp.replace(path)


def align_rows(
    input_csv: str | Path,
    output_dir: str | Path,
    *,
    usalign_path: str | Path = DEFAULT_USALIGN,
    workers: int = 1,
    timeout: int = 300,
    keep_raw: bool = False,
    limit: int | None = None,
) -> dict[str, Path]:
    input_path = Path(input_csv).expanduser().resolve()
    out = Path(output_dir).expanduser().resolve()
    out.mkdir(parents=True, exist_ok=True)
    rows = _read_input_csv(input_path)
    if limit is not None:
        rows = rows[:limit]
    usalign = Path(usalign_path).expanduser().resolve()
    results: list[tuple[dict[str, object], list[dict[str, object]], list[dict[str, object]]]] = []
    if workers <= 1:
        for row in rows:
            results.append(align_one_row(row, input_base=input_path.parent, output_dir=out, usalign_path=usalign, timeout=timeout, keep_raw=keep_raw))
    else:
        with ThreadPoolExecutor(max_workers=workers) as executor:
            futures = [executor.submit(align_one_row, row, input_base=input_path.parent, output_dir=out, usalign_path=usalign, timeout=timeout, keep_raw=keep_raw) for row in rows]
            for future in as_completed(futures):
                results.append(future.result())
        order = {row["design_id"]: index for index, row in enumerate(rows)}
        results.sort(key=lambda item: order.get(str(item[0].get("design_id", "")), len(order)))
    summaries = [item[0] for item in results]
    chain_rows = [row for item in results for row in item[1]]
    motif_rows = [row for item in results for row in item[2]]
    summary_path = out / "alignment_summary.csv"
    chain_path = out / "chain_pairs.csv"
    motif_path = out / "motifs.csv"
    _write_rows(summary_path, summaries, preferred=("design_id", "sequence_id", "status", "error"))
    _write_rows(chain_path, chain_rows, preferred=("design_id", "backbone_chain", "predicted_chain", "status"))
    _write_rows(motif_path, motif_rows, preferred=("design_id", "motif_name", "status", "error"))
    return {"summary": summary_path, "chain_pairs": chain_path, "motifs": motif_path}


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Generic v37 backbone/design structure alignment")
    parser.add_argument("--input-csv", required=True, help="CSV with design_id, backbone_structure_path, predicted_structure_path")
    parser.add_argument("--outdir", required=True, help="output directory for three CSV tables")
    parser.add_argument("--usalign", default=DEFAULT_USALIGN, help="USalign executable")
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--timeout", type=int, default=300, help="timeout per USalign call in seconds")
    parser.add_argument("--keep-raw", action="store_true", help="retain per-design raw USalign text under outdir/raw")
    parser.add_argument("--limit", type=int, help="process only the first N input rows")
    args = parser.parse_args(argv)
    if args.workers < 1:
        parser.error("--workers must be >= 1")
    paths = align_rows(
        args.input_csv,
        args.outdir,
        usalign_path=args.usalign,
        workers=args.workers,
        timeout=args.timeout,
        keep_raw=args.keep_raw,
        limit=args.limit,
    )
    print("alignment complete")
    for name, path in paths.items():
        print(f"{name}: {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
