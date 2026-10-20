#!/usr/bin/env python3
"""Model-neutral confidence metric extraction from full-data payloads.

Adapters provide the raw JSON/NPZ path; this module only performs deterministic
CPU reductions.  Missing or incompatible arrays result in omitted metrics,
never fabricated zeros.
"""

from __future__ import annotations

import math
from typing import Any, Iterable

from .ipsae import pairwise_ipae_ipsae


def _numbers(value: Any) -> list[float]:
    if hasattr(value, "tolist"):
        value = value.tolist()
    if isinstance(value, (int, float)) and math.isfinite(float(value)):
        return [float(value)]
    if isinstance(value, (list, tuple)):
        output: list[float] = []
        for item in value:
            output.extend(_numbers(item))
        return output
    return []


def _matrix(value: Any) -> list[list[float]]:
    if hasattr(value, "tolist"):
        value = value.tolist()
    if not isinstance(value, (list, tuple)) or not value:
        return []
    result: list[list[float]] = []
    for row in value:
        values = _numbers(row)
        if values:
            result.append(values)
    width = len(result[0]) if result else 0
    return [row[:width] for row in result if len(row) == width] if width else []


def _mean(values: Iterable[float]) -> float | str:
    values = [value for value in values if math.isfinite(value)]
    return sum(values) / len(values) if values else ""


def _pick(data: dict[str, Any], names: tuple[str, ...]) -> Any:
    lowered = {str(key).lower(): value for key, value in data.items()}
    for name in names:
        if name.lower() in lowered:
            return lowered[name.lower()]
    return None


def extract_confidence_metrics(data: dict[str, Any]) -> dict[str, object]:
    """Return scalar, per-chain and per-chain-pair confidence reductions."""
    output: dict[str, object] = {}
    plddt = _numbers(_pick(data, ("plddt", "atom_plddts", "token_plddts")))
    if plddt:
        output["mean_plddt"] = _mean(plddt)
    pae = _matrix(_pick(data, ("token_pair_pae", "pae", "predicted_aligned_error")))
    chain_ids_raw = _pick(data, ("token_asym_id", "token_chain_ids", "chain_id", "asym_id"))
    if hasattr(chain_ids_raw, "tolist"):
        chain_ids_raw = chain_ids_raw.tolist()
    chain_ids = [str(x) for x in chain_ids_raw] if isinstance(chain_ids_raw, (list, tuple)) else []
    if pae:
        values = [value for row in pae for value in row]
        output["mean_pae"] = _mean(values)
        if chain_ids and len(chain_ids) == len(pae):
            chain_metrics: dict[str, dict[str, object]] = {}
            for chain in sorted(set(chain_ids)):
                indices = [i for i, value in enumerate(chain_ids) if value == chain]
                chain_metrics[chain] = {
                    "mean_pae": _mean(pae[i][j] for i in indices for j in indices),
                    "mean_plddt": _mean(plddt[i] for i in indices if i < len(plddt)),
                    "length": len(indices),
                }
            pair_scores = pairwise_ipae_ipsae(pae, chain_ids)
            pair_metrics: dict[str, dict[str, object]] = {
                pair: {"mean_interface_pae": values[0], "ipsae": values[1]}
                for pair, values in pair_scores.items()
            }
            output["chain_metrics"] = chain_metrics
            output["chain_pair_metrics"] = pair_metrics
            cross = [pae[i][j] for i, a in enumerate(chain_ids) for j, b in enumerate(chain_ids) if a != b]
            if cross:
                output["mean_interface_pae"] = _mean(cross)
            if pair_scores:
                output["ipsae"] = sum(values[1] for values in pair_scores.values()) / len(pair_scores)
    # IPSAE is model-specific.  Preserve a supplied value and do not derive a
    # different score from PAE without the model's documented normalization.
    ipsae = _pick(data, ("ipsae", "overall_ipsae_mean", "interface_psae"))
    if ipsae not in (None, ""):
        output["ipsae"] = ipsae
    return output
