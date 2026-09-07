#!/usr/bin/env python3
"""Build run-local model inputs without changing user-provided source files."""

from __future__ import annotations

import copy
import json
import shutil
from pathlib import Path
from typing import Any, Mapping

from .manifest import metadata, read_fasta, sequence_from_row


def copy_or_none(source: str | object, destination: str | Path) -> str:
    if not source:
        return ""
    source_path = Path(str(source)).expanduser().resolve()
    if not source_path.is_file():
        raise FileNotFoundError(source_path)
    output = Path(destination)
    output.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(source_path, output)
    return str(output.resolve())


def protein_records(row: Mapping[str, object]) -> list[tuple[str, str]]:
    """Return ordered ``(chain_id, sequence)`` pairs for model input builders.

    The compact v37 manifest accepts either one sequence, ``A:B``/``A|B``
    separated chains, or a FASTA.  ``chain_ids`` is optional; when omitted,
    chains are named A, B, ... in input order.  Complex inputs that contain a
    ligand or non-protein entity must be supplied explicitly through the
    model-specific input path fields.
    """

    sequence = sequence_from_row(row).strip().upper()
    if sequence:
        chains = [part.strip() for part in sequence.replace("|", ":").split(":") if part.strip()]
        raw_ids = str(row.get("chain_ids", "")).replace(";", ",").split(",")
        chain_ids = [item.strip() for item in raw_ids if item.strip()]
        if chain_ids and len(chain_ids) != len(chains):
            raise ValueError(
                f"design {row.get('design_id', '')}: chain_ids count does not match sequence chains"
            )
        if not chain_ids:
            chain_ids = [chr(ord("A") + index) for index in range(len(chains))]
        if len(set(chain_ids)) != len(chain_ids):
            raise ValueError(f"design {row.get('design_id', '')}: duplicate chain_ids")
        return list(zip(chain_ids, chains))

    fasta_path = str(row.get("fasta_path", "")).strip()
    if fasta_path:
        records = read_fasta(fasta_path)
        if not records:
            raise ValueError(f"FASTA contains no sequence: {fasta_path}")
        output: list[tuple[str, str]] = []
        for index, (name, value) in enumerate(records):
            chain_id = name.rsplit("__chain__", 1)[-1] if "__chain__" in name else chr(ord("A") + index)
            output.append((chain_id, value))
        return output
    raise ValueError(f"design {row.get('design_id', '')}: sequence/fasta_path is required")


def _protein_entries(payload: Mapping[str, Any]) -> list[dict[str, Any]]:
    entries: list[dict[str, Any]] = []
    for item in payload.get("sequences", []) if isinstance(payload.get("sequences"), list) else []:
        if isinstance(item, dict) and isinstance(item.get("protein"), dict):
            entries.append(item["protein"])
    return entries


def _protein_chain_ids(proteins: list[Mapping[str, Any]]) -> list[list[str]]:
    """Return template chain IDs while rejecting ambiguous AF3 templates."""

    output: list[list[str]] = []
    seen: set[str] = set()
    for index, protein in enumerate(proteins):
        raw_ids = protein.get("id")
        if isinstance(raw_ids, list):
            chain_ids = [str(item).strip() for item in raw_ids if str(item).strip()]
        elif raw_ids not in (None, ""):
            chain_ids = [str(raw_ids).strip()]
        else:
            chain_ids = []
        if not chain_ids:
            raise ValueError(f"AF3 protein entry {index} has no chain id")
        if len(set(chain_ids)) != len(chain_ids):
            raise ValueError(f"AF3 protein entry {index} has duplicate chain ids")
        overlap = sorted(set(chain_ids) & seen)
        if overlap:
            raise ValueError(f"AF3 template reuses chain id(s): {', '.join(overlap)}")
        seen.update(chain_ids)
        output.append(chain_ids)
    return output


def _sequence(value: object, label: str) -> str:
    sequence = str(value).strip().replace(" ", "").replace("\n", "").upper()
    if not sequence:
        raise ValueError(f"{label} is empty")
    if ":" in sequence or "|" in sequence:
        raise ValueError(f"{label} must contain one protein-chain sequence")
    return sequence


def _af3_sequence_map(
    row: Mapping[str, object],
    proteins: list[dict[str, Any]],
    sequence_spec: object,
) -> dict[str, str]:
    """Resolve every AF3 protein chain to a replacement sequence."""

    chain_groups = _protein_chain_ids(proteins)
    template_ids = [chain_id for group in chain_groups for chain_id in group]
    expected = set(template_ids)

    if isinstance(sequence_spec, dict):
        provided = {str(key).strip(): value for key, value in sequence_spec.items()}
        provided_ids = {key for key in provided if key}
        missing = sorted(expected - provided_ids)
        extra = sorted(provided_ids - expected)
        if missing or extra:
            details = []
            if missing:
                details.append(f"missing chain id(s): {', '.join(missing)}")
            if extra:
                details.append(f"unknown chain id(s): {', '.join(extra)}")
            raise ValueError("AF3 sequence mapping does not match template; " + "; ".join(details))
        result = {chain_id: _sequence(provided[chain_id], f"AF3 sequence for chain {chain_id}") for chain_id in expected}
    elif isinstance(sequence_spec, list):
        if len(sequence_spec) != len(proteins):
            raise ValueError(
                "AF3 sequence list count does not match protein entry count "
                f"({len(sequence_spec)} != {len(proteins)})"
            )
        result = {}
        for chain_ids, value in zip(chain_groups, sequence_spec):
            sequence = _sequence(value, f"AF3 sequence for chain(s) {','.join(chain_ids)}")
            result.update({chain_id: sequence for chain_id in chain_ids})
    else:
        records = protein_records(row)
        provided = {chain_id: sequence for chain_id, sequence in records}
        provided_ids = set(provided)
        missing = sorted(expected - provided_ids)
        extra = sorted(provided_ids - expected)
        if missing or extra:
            details = []
            if missing:
                details.append(f"missing chain id(s): {', '.join(missing)}")
            if extra:
                details.append(f"unknown chain id(s): {', '.join(extra)}")
            raise ValueError("row sequence chains do not match AF3 template; " + "; ".join(details))
        result = {chain_id: _sequence(provided[chain_id], f"row sequence for chain {chain_id}") for chain_id in expected}

    for chain_ids in chain_groups:
        sequences = {result[chain_id] for chain_id in chain_ids}
        if len(sequences) != 1:
            raise ValueError(
                "AF3 protein entry has multiple chain IDs with different sequences: "
                + ", ".join(chain_ids)
            )
    return result


_NO_MSA_KEYS = {
    "msa",
    "unpairedmsa",
    "pairedmsa",
    "templatemsa",
    "templates",
    "a3m",
    "a3mpath",
    "msapath",
    "msafile",
    "msafilepath",
}
_EMPTY_MSA_STRINGS = {"", "empty", "none", "null"}


def _is_empty_msa(value: object) -> bool:
    if value is None or value is False:
        return True
    if isinstance(value, str):
        return value.strip().lower() in _EMPTY_MSA_STRINGS
    if isinstance(value, (list, tuple, dict, set)):
        return len(value) == 0
    return False


def _find_nonempty_msa(value: object, path: str = "input") -> str | None:
    if isinstance(value, Mapping):
        for key, child in value.items():
            key_text = str(key)
            normalized = key_text.replace("_", "").replace("-", "").lower()
            child_path = f"{path}.{key_text}"
            if normalized in _NO_MSA_KEYS and not _is_empty_msa(child):
                return child_path
            violation = _find_nonempty_msa(child, child_path)
            if violation:
                return violation
    elif isinstance(value, list):
        for index, child in enumerate(value):
            violation = _find_nonempty_msa(child, f"{path}[{index}]")
            if violation:
                return violation
    return None


def _validate_no_msa_input(source: str | Path, model: str) -> None:
    source_path = Path(str(source)).expanduser().resolve()
    if not source_path.is_file():
        raise FileNotFoundError(source_path)
    try:
        if model == "protenix":
            payload = json.loads(source_path.read_text(encoding="utf-8"))
        elif model == "boltz2":
            import yaml

            payload = yaml.safe_load(source_path.read_text(encoding="utf-8"))
        else:
            raise ValueError(f"unsupported no-MSA input model: {model}")
    except (OSError, json.JSONDecodeError) as exc:
        raise ValueError(f"cannot parse {model} input: {source_path}") from exc
    except Exception as exc:
        if model == "boltz2":
            raise ValueError(f"cannot parse Boltz-2 YAML input: {source_path}") from exc
        raise
    if not isinstance(payload, (Mapping, list)):
        raise ValueError(f"{model} input must contain an object or list: {source_path}")
    violation = _find_nonempty_msa(payload)
    if violation:
        raise ValueError(f"{model} input contains non-empty MSA at {violation}: {source_path}")


def build_af3_json(row: Mapping[str, object], destination: str | Path, *, seed: int) -> str:
    """Copy a ready AF3 JSON or fill a simple AF3 protein template.

    A row may provide ``af3_json_path`` directly.  For a template in
    ``target_spec_path``, ``metadata_json.af3_sequences`` can be either
    ``{"A": "SEQUENCE"}`` or a list of sequences.  If the template has one
    protein entity, the manifest ``sequence`` is used automatically.
    """

    source = row.get("af3_json_path") or row.get("target_spec_path")
    if not source:
        raise FileNotFoundError("af3_json_path/target_spec_path is not provided")
    source_path = Path(str(source)).expanduser().resolve()
    payload = json.loads(source_path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise ValueError(f"AF3 JSON must contain an object: {source_path}")
    payload = copy.deepcopy(payload)
    design_id = str(row["design_id"])
    payload["name"] = design_id
    payload["modelSeeds"] = [int(seed)]
    proteins = _protein_entries(payload)
    meta = metadata(row)
    if not proteins:
        raise ValueError(f"AF3 JSON contains no protein entities: {source_path}")
    sequence_map = _af3_sequence_map(row, proteins, meta.get("af3_sequences"))
    for protein, chain_ids in zip(proteins, _protein_chain_ids(proteins)):
        protein["sequence"] = sequence_map[chain_ids[0]]
    for protein in proteins:
        for key in ("unpairedMsa", "pairedMsa", "templateMsa", "templates"):
            if key in protein and meta.get("clear_msa", True):
                protein[key] = []
    output = Path(destination)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    return str(output.resolve())


def _copy_or_build_model_input(
    row: Mapping[str, object],
    destination: str | Path,
    field: str,
    builder: Any,
    model: str | None = None,
) -> str:
    source = row.get(field)
    if source:
        if model:
            _validate_no_msa_input(str(source), model)
        return copy_or_none(source, destination)
    output = Path(destination)
    output.parent.mkdir(parents=True, exist_ok=True)
    builder(output)
    return str(output.resolve())


def build_protenix_json(row: Mapping[str, object], destination: str | Path) -> str:
    """Build the no-MSA Protenix JSON accepted by the local wrapper."""

    records = protein_records(row)

    def builder(output: Path) -> None:
        payload = [
            {
                "name": str(row["design_id"]),
                "sequences": [
                    {
                        "proteinChain": {
                            "count": 1,
                            "id": [chain_id],
                            "sequence": sequence,
                        }
                    }
                    for chain_id, sequence in records
                ],
            }
        ]
        output.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    return _copy_or_build_model_input(row, destination, "protenix_input_path", builder, "protenix")


def build_boltz_yaml(row: Mapping[str, object], destination: str | Path) -> str:
    """Build a minimal Boltz no-MSA YAML, unless an explicit input is given."""

    records = protein_records(row)

    def builder(output: Path) -> None:
        lines = ["version: 1", "sequences:"]
        for chain_id, sequence in records:
            lines.extend(
                [
                    "  - protein:",
                    f"      id: [{chain_id}]",
                    f"      sequence: {sequence}",
                    "      msa: empty",
                ]
            )
        output.write_text("\n".join(lines) + "\n", encoding="utf-8")

    return _copy_or_build_model_input(row, destination, "boltz_input_path", builder, "boltz2")


def build_opendde_fasta(row: Mapping[str, object], destination: str | Path) -> str:
    """Build an OpenDDE FASTA (colon-separated chains), unless supplied."""

    records = protein_records(row)

    def builder(output: Path) -> None:
        output.write_text(
            f">{row['design_id']}\n{':'.join(sequence for _chain, sequence in records)}\n",
            encoding="utf-8",
        )

    return _copy_or_build_model_input(row, destination, "opendde_input_path", builder)


def status_from_models(statuses: Mapping[str, str]) -> str:
    values = list(statuses.values())
    if not values:
        return "missing"
    if all(value == "success" for value in values):
        return "success"
    if any(value in {"success", "partial"} for value in values):
        return "partial"
    if all(value == "missing" for value in values):
        return "missing"
    return "failed"
