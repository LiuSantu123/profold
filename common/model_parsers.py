#!/usr/bin/env python3
"""Small, defensive parsers for model output contracts.

These parsers intentionally do not import the historical monitoring scripts.
They preserve raw paths and summary values while treating missing detail files
as ``partial`` instead of manufacturing confidence values.
"""

from __future__ import annotations

import csv
import hashlib
import json
import math
from pathlib import Path
from typing import Any, Iterable


def _number(value: object) -> object:
    if value in (None, "", "N/A", "NA", "null"):
        return ""
    try:
        number = float(value)  # type: ignore[arg-type]
        if not math.isfinite(number):
            return ""
        return int(number) if number.is_integer() else number
    except (TypeError, ValueError):
        return value


def _scaled_plddt(value: object) -> object:
    """Return pLDDT on the common 0--100 scale.

    Protenix/AF3 summaries normally use 0--100, while Boltz-2 and some
    OpenDDE fields use 0--1.  The threshold mirrors the monitor memory and
    deliberately leaves values above 1.1 unchanged.
    """

    number = _number(value)
    if isinstance(number, (int, float)) and 0 <= number <= 1.1:
        return _number(number * 100)
    return number


def _files(root: Path, patterns: Iterable[str]) -> list[Path]:
    found: list[Path] = []
    seen: set[Path] = set()
    for pattern in patterns:
        for path in root.rglob(pattern):
            if path.is_file() and path not in seen:
                found.append(path)
                seen.add(path)
    return sorted(found)


def _choose(paths: list[Path], design_id: str) -> Path | None:
    if not paths:
        return None

    def is_exact(path: Path) -> bool:
        if design_id in path.parts:
            return True
        name = path.name
        if name == design_id or path.stem == design_id:
            return True
        return any(
            name.startswith(f"{design_id}{separator}")
            for separator in (".", "_", "-", "__")
        )

    exact = [path for path in paths if is_exact(path)]
    if exact:
        return sorted(
            exact,
            key=lambda item: (0 if item.name.startswith(f"{design_id}_") else 1, str(item)),
        )[0]
    # A single-output directory is safe to associate with the sole input.  In
    # a multi-design directory, refusing an ambiguous match is safer than
    # attaching another design's structure or confidence score.
    return paths[0] if len(paths) == 1 else None


def _load_json(path: Path | None) -> dict[str, Any]:
    if path is None:
        return {}
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return {}
    return value if isinstance(value, dict) else {}


def _structure(root: Path, design_id: str) -> Path | None:
    paths = _files(root, (f"{design_id}*.cif", f"{design_id}*.mmcif", f"{design_id}*.pdb", "*.cif", "*.mmcif", "*.pdb"))
    return _choose(paths, design_id)


def _status(summary: Path | None, structure: Path | None, root: Path) -> tuple[str, str]:
    if summary and structure:
        return "success", ""
    if summary or structure:
        return "partial", "summary or structure artifact is missing"
    if root.exists():
        return "failed", "no recognized model output found"
    return "missing", "model output directory does not exist"


def _summary_row(path: Path | None, *, name_keys: tuple[str, ...] = ("name", "design_id", "model_name")) -> dict[str, Any]:
    if path is None:
        return {}
    data = _load_json(path)
    return data


def _csv_rows(path: Path | None) -> list[dict[str, str]]:
    if path is None:
        return []
    try:
        with path.open(encoding="utf-8", newline="") as handle:
            return [dict(row) for row in csv.DictReader(handle)]
    except (OSError, csv.Error):
        return []


def _csv_row(path: Path | None, design_id: str) -> dict[str, str]:
    rows = _csv_rows(path)
    if not rows:
        return {}
    id_fields = ("design_id", "sid", "id", "name", "sequence_id")
    for row in rows:
        for field in id_fields:
            if row.get(field, "").strip() == design_id:
                return row
    return rows[0] if len(rows) == 1 else {}


def _value(row: dict[str, Any], names: Iterable[str]) -> object:
    lowered = {str(key).lower(): value for key, value in row.items()}
    for name in names:
        if name.lower() in lowered:
            return lowered[name.lower()]
    return ""


def parse_esmfold(
    root: str | Path,
    design_id: str,
    *,
    prefix: str = "esmfold",
    native: bool = False,
) -> dict[str, object]:
    """Parse native ESMFold PDB or a custom predictor's summary/structure.

    Level1 calls this with the historical ``esmfold`` prefix.  ``prefix`` is
    retained as a small compatibility hook for callers that need another
    output column namespace.
    """

    directory = Path(root).resolve()
    if native:
        # Exact native filename prevents cross-design and legacy CIF reuse.
        pdb = directory / f"{design_id}.pdb"
        result = {
            f"{prefix}_status": "failed" if directory.exists() else "missing",
            f"{prefix}_error": "native ESMFold PDB is missing",
            f"{prefix}_output_dir": str(directory) if directory.exists() else "",
            f"{prefix}_summary_path": "",
            f"{prefix}_structure_path": "",
        }
        if pdb.is_file():
            try:
                import gemmi

                structure = gemmi.read_structure(str(pdb))
                atoms = [atom for model in structure for chain in model for residue in chain for atom in residue]
                if not atoms or not any(atom.name == "CA" for atom in atoms):
                    raise ValueError("PDB has no protein coordinates")
                scores = [float(atom.b_iso) for atom in atoms]
                if any(not math.isfinite(score) or not 0 <= score <= 100 for score in scores):
                    raise ValueError("invalid PDB pLDDT")
                if any(not math.isfinite(v) for atom in atoms for v in (atom.pos.x, atom.pos.y, atom.pos.z)):
                    raise ValueError("invalid PDB coordinates")
                result.update({
                    f"{prefix}_status": "success",
                    f"{prefix}_error": "",
                    f"{prefix}_structure_path": str(pdb),
                    f"{prefix}_mean_plddt": sum(scores) / len(scores),
                })
            except Exception as exc:
                result[f"{prefix}_status"] = "failed"
                result[f"{prefix}_error"] = f"invalid native ESMFold PDB: {exc}"
        return result
    summary_path = _choose(_files(directory, ("esmfold2_summary.csv", "*_summary.csv", "*.csv")), design_id)
    metrics: dict[str, object] = {}
    if summary_path:
        source = _csv_row(summary_path, design_id)
        for source_key, target_key in (
            ("mean_plddt", f"{prefix}_mean_plddt"),
            ("plddt", f"{prefix}_mean_plddt"),
            ("ptm", f"{prefix}_ptm"),
            ("iptm", f"{prefix}_iptm"),
        ):
            if source_key in source and target_key not in metrics:
                metrics[target_key] = _number(source[source_key])
    structure = _structure(directory, design_id)
    status, error = _status(summary_path, structure, directory)
    result: dict[str, object] = {
        f"{prefix}_status": status,
        f"{prefix}_error": error,
        f"{prefix}_output_dir": str(directory) if directory.exists() else "",
        f"{prefix}_summary_path": str(summary_path) if summary_path else "",
        f"{prefix}_structure_path": str(structure) if structure else "",
    }
    result.update(metrics)
    return result


def _af3_detail_path(summary: Path | None) -> Path | None:
    if summary is None:
        return None
    candidate = summary.with_name(
        summary.name.replace("_summary_confidences", "_confidences").replace(
            "_summary_confidence", "_confidence"
        )
    )
    return candidate if candidate.is_file() else None


def _af3_variant_match(path: Path, design_id: str, variant: str) -> bool:
    base = f"{design_id}_af3_{variant}"
    return path.name == base or path.name.startswith(f"{base}_") or path.name.startswith(f"{base}.")


def _af3_variant_files(directory: Path, design_id: str, variant: str, patterns: Iterable[str]) -> list[Path]:
    return [
        path
        for path in _files(directory, patterns)
        if _af3_variant_match(path, design_id, variant)
    ]


def _mean_value(value: object) -> object:
    if isinstance(value, (list, tuple)):
        numbers: list[float] = []
        for item in value:
            parsed = _number(item)
            if isinstance(parsed, (int, float)):
                numbers.append(float(parsed))
        return _number(sum(numbers) / len(numbers)) if numbers else ""
    return _number(value)


def _parse_af3_variant(
    directory: Path,
    design_id: str,
    prefix: str,
    *,
    variant: str | None,
) -> dict[str, object]:
    summary_patterns = ("*_summary_confidences.json", "*_summary_confidence*.json", "confidence.json")
    structure_patterns = (
        f"{design_id}*.cif",
        f"{design_id}*.mmcif",
        f"{design_id}*.pdb",
        "*.cif",
        "*.mmcif",
        "*.pdb",
    )
    if variant is None:
        summary = _choose(_files(directory, summary_patterns), design_id)
        structure = _structure(directory, design_id)
    else:
        summary = _choose(
            _af3_variant_files(directory, design_id, variant, summary_patterns),
            design_id,
        )
        structure = _choose(
            _af3_variant_files(directory, design_id, variant, structure_patterns),
            design_id,
        )
    detail = _af3_detail_path(summary)
    data = _summary_row(summary)
    detail_data = _load_json(detail)
    status, error = _status(summary, structure, directory)
    result: dict[str, object] = {
        f"{prefix}_status": status,
        f"{prefix}_error": error,
        f"{prefix}_output_dir": str(directory) if directory.exists() else "",
        f"{prefix}_summary_json": str(summary) if summary else "",
        f"{prefix}_detail_json": str(detail) if detail else "",
        f"{prefix}_structure_path": str(structure) if structure else "",
    }
    aliases = {
        "iptm": "iptm",
        "ptm": "ptm",
        "ranking_score": "ranking_score",
        "fraction_disordered": "fraction_disordered",
        "has_clash": "has_clash",
    }
    for source, suffix in aliases.items():
        if source in data:
            result[f"{prefix}_{suffix}"] = _number(data[source])

    for source in ("plddt", "atom_plddts", "mean_plddt"):
        if source in detail_data:
            result[f"{prefix}_mean_plddt"] = _mean_value(detail_data[source])
            break
    if f"{prefix}_mean_plddt" not in result:
        for source in ("mean_plddt", "plddt"):
            if source in data:
                result[f"{prefix}_mean_plddt"] = _mean_value(data[source])
                break
    if variant is not None:
        model_name = f"{design_id}_af3_{variant}"
        for record_path in directory.rglob("result.json"):
            if record_path.parent.name != model_name:
                continue
            record = _load_json(record_path)
            try:
                if not record.get("artifacts"):
                    raise ValueError("archive has no artifact manifest")
                for item in record["artifacts"]:
                    artifact = record_path.parent / item["name"]
                    if artifact.parent != record_path.parent or hashlib.sha256(artifact.read_bytes()).hexdigest() != item["sha256"]:
                        raise ValueError("archive artifact checksum mismatch")
            except (OSError, ValueError, KeyError, TypeError) as exc:
                result[f"{prefix}_status"] = "failed"
                result[f"{prefix}_error"] = str(exc)
                return result
            archived = record.get("result", {})
            if archived.get("model_name") == model_name and archived.get("status") == "success":
                for key, value in archived.items():
                    if key not in {"status", "error", "model_name", "af3_cif_path"}:
                        result[f"{prefix}_{key}"] = value if value is not None else ""
    return result


def _combined_status(statuses: Iterable[str]) -> str:
    values = list(statuses)
    if values and all(value == "success" for value in values):
        return "success"
    if any(value in {"success", "partial"} for value in values):
        return "partial"
    if values and all(value == "missing" for value in values):
        return "missing"
    return "failed"


def parse_af3(root: str | Path, design_id: str, *, mode: str = "direct") -> dict[str, object]:
    """Parse direct AF3 output or the two outputs produced by ``for_wj``.

    ``for_wj_cid`` archives the with-ligand and without-ligand artifacts in
    separate flat directories after its monitor cleans the task folders.  The
    variant-aware lookup therefore matches the full output basename instead of
    choosing one arbitrary JSON/CIF pair.
    """

    directory = Path(root).resolve()
    if mode == "direct":
        return _parse_af3_variant(directory, design_id, "af3", variant=None)
    if mode not in {"for_wj_cid", "cid"}:
        raise ValueError(f"unsupported AF3 mode: {mode}")

    variants = {
        variant: _parse_af3_variant(directory, design_id, f"af3_{variant}", variant=variant)
        for variant in ("with_lig", "without_lig")
    }
    result: dict[str, object] = {
        "af3_output_dir": str(directory) if directory.exists() else "",
        "af3_variants": "with_lig,without_lig",
    }
    for variant_result in variants.values():
        result.update(variant_result)

    statuses = [str(variants[variant][f"af3_{variant}_status"]) for variant in variants]
    result["af3_status"] = _combined_status(statuses)
    errors = [
        f"{variant}: {variants[variant][f'af3_{variant}_error']}"
        for variant in variants
        if variants[variant][f"af3_{variant}_error"]
    ]
    result["af3_error"] = "; ".join(errors)

    # Keep the old af3_* contract usable by selectors written for direct AF3.
    # with_lig is the primary result for CID screening; fall back to the other
    # variant when it is the only artifact that exists.
    primary_variant = next(
        (
            variant
            for variant in ("with_lig", "without_lig")
            if variants[variant][f"af3_{variant}_summary_json"]
            or variants[variant][f"af3_{variant}_structure_path"]
        ),
        "with_lig",
    )
    primary = variants[primary_variant]
    if mode == "cid":
        primary_variant = "with_lig"
        primary = variants[primary_variant]
    for key, value in primary.items():
        marker = f"af3_{primary_variant}_"
        if key.startswith(marker):
            result[f"af3_{key[len(marker):]}"] = value
    result["af3_status"] = _combined_status(statuses)
    result["af3_error"] = "; ".join(errors)
    for state, variant in (("apo", "without_lig"), ("holo", "with_lig")):
        marker = f"af3_{variant}_"
        for key, value in variants[variant].items():
            if key.startswith(marker):
                result[f"af3_{state}_{key[len(marker):]}"] = value
    result["af3_states"] = "apo,holo"
    return result


def parse_summary_model(root: str | Path, design_id: str, model: str) -> dict[str, object]:
    directory = Path(root).resolve()
    if model == "protenix":
        patterns = ("*_summary_confidence_sample_*.json", "*_summary_confidences.json")
    elif model == "boltz2":
        patterns = ("confidence_*.json",)
    elif model == "opendde":
        patterns = ("*_summary_confidence_sample_*.json", "*_summary_confidences.json")
    else:
        raise ValueError(f"unsupported summary model: {model}")
    summary = _choose(_files(directory, patterns), design_id)
    data = _summary_row(summary)
    structure = _structure(directory, design_id)
    status, error = _status(summary, structure, directory)
    prefix = model
    result: dict[str, object] = {
        f"{prefix}_status": status,
        f"{prefix}_error": error,
        f"{prefix}_output_dir": str(directory) if directory.exists() else "",
        f"{prefix}_summary_json": str(summary) if summary else "",
        f"{prefix}_structure_path": str(structure) if structure else "",
    }
    if model in {"protenix", "opendde"}:
        full_patterns = ("*_full_data_sample_*.json", "*_full_data*.json")
        full_data = _choose(_files(directory, full_patterns), design_id)
        result[f"{prefix}_full_data_json"] = str(full_data) if full_data else ""
    if model == "boltz2":
        affinity = _choose(_files(directory, ("affinity_*.json",)), design_id)
        affinity_data = _load_json(affinity)
        result[f"{prefix}_affinity_json"] = str(affinity) if affinity else ""
        for key in ("affinity_pred_value", "affinity_probability_binary"):
            if key in affinity_data:
                result[f"{prefix}_{key}"] = _number(affinity_data[key])
    for key in (
        "iptm",
        "ptm",
        "ranking_score",
        "plddt",
        "confidence_score",
        "complex_plddt",
        "complex_iplddt",
        "complex_pde",
        "complex_ipde",
        "ligand_iptm",
        "protein_iptm",
        "affinity_pred_value",
        "affinity_probability_binary",
        "has_clash",
        "disorder",
        "gpde",
        "chain_ptm",
        "chain_iptm",
        "chain_pair_iptm",
        "chain_plddt",
        "chain_gpde",
        "chain_pair_gpde",
    ):
        if key in data:
            result[f"{prefix}_{key}"] = _number(data[key])
    if "plddt" in data:
        result[f"{prefix}_mean_plddt"] = _scaled_plddt(data["plddt"])
    elif "complex_plddt" in data:
        result[f"{prefix}_mean_plddt"] = _scaled_plddt(data["complex_plddt"])
    if model == "boltz2" and "complex_iplddt" in data:
        result[f"{prefix}_mean_interface_plddt"] = _scaled_plddt(data["complex_iplddt"])
    return result


def parse_netsolp(path: str | Path, design_id: str) -> dict[str, object]:
    """Parse the extracted NetSolP CSV produced by ``run_netsolp.sh``."""

    input_path = Path(path).resolve()
    csv_path = input_path
    if input_path.is_dir():
        candidates = _files(input_path, ("*_extracted.csv", "*.csv"))
        csv_path = _choose(candidates, design_id) or (candidates[0] if len(candidates) == 1 else input_path / "")
    row = _csv_row(csv_path if csv_path.is_file() else None, design_id)
    solubility = _number(_value(row, ("predicted_solubility", "solubility", "avg_solubility")))
    usability = _number(_value(row, ("predicted_usability", "usability", "avg_usability")))
    has_metric = solubility != "" or usability != ""
    if row and has_metric:
        status, error = "success", ""
    elif row:
        status, error = "partial", "design row found but no solubility metric"
    elif csv_path.is_file():
        status, error = "failed", "design_id not found in NetSolP CSV"
    elif input_path.exists():
        status, error = "failed", "NetSolP CSV is missing"
    else:
        status, error = "missing", "NetSolP output does not exist"
    return {
        "netsolp_status": status,
        "netsolp_error": error,
        "netsolp_csv_path": str(csv_path) if csv_path.is_file() else "",
        "netsolp_predicted_solubility": solubility,
        "netsolp_predicted_usability": usability,
    }


def parse_temberture(path: str | Path, design_id: str) -> dict[str, object]:
    """Parse the ``tem_TM``/``tem_SC`` CSV from the historical TemBERTure script."""

    input_path = Path(path).resolve()
    csv_path = input_path
    if input_path.is_dir():
        candidates = _files(input_path, ("*temberture*.csv", "*.csv"))
        csv_path = _choose(candidates, design_id) or (candidates[0] if len(candidates) == 1 else input_path / "")
    row = _csv_row(csv_path if csv_path.is_file() else None, design_id)
    tm = _number(_value(row, ("tem_TM", "tm", "predicted_tm", "temperature")))
    score = _number(_value(row, ("tem_SC", "thermostability_score", "thermophilicity_score", "score")))
    has_metric = tm != "" or score != ""
    if row and has_metric:
        status, error = "success", ""
    elif row:
        status, error = "partial", "design row found but no thermal metric"
    elif csv_path.is_file():
        status, error = "failed", "design_id not found in TemBERTure CSV"
    elif input_path.exists():
        status, error = "failed", "TemBERTure CSV is missing"
    else:
        status, error = "missing", "TemBERTure output does not exist"
    return {
        "temberture_status": status,
        "temberture_error": error,
        "temberture_csv_path": str(csv_path) if csv_path.is_file() else "",
        "temberture_tm": tm,
        "temberture_score": score,
    }


def unavailable(model: str, reason: str) -> dict[str, object]:
    return {f"{model}_status": "missing", f"{model}_error": reason}
