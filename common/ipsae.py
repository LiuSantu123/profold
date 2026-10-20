#!/usr/bin/env python3
"""Small, dependency-free implementation of the DunbrackLab ipSAE reduction."""

from __future__ import annotations

from math import isfinite


def d0(length: int) -> float:
    return max(1.0, 1.24 * (length - 15) ** (1.0 / 3.0) - 1.8) if length > 27 else 1.0


def pairwise_ipae_ipsae(pae: list[list[float]], chain_ids: list[str], cutoff: float = 15.0) -> dict[str, tuple[float, float]]:
    """Return mean inter-chain PAE and directional-max ipSAE per chain pair."""
    unique = list(dict.fromkeys(str(x) for x in chain_ids))
    result: dict[str, tuple[float, float]] = {}
    for ai, left in enumerate(unique):
        li = [i for i, value in enumerate(chain_ids) if str(value) == left]
        for right in unique[ai + 1:]:
            ri = [i for i, value in enumerate(chain_ids) if str(value) == right]
            scale = d0(len(li) + len(ri))
            blocks: list[float] = []
            directional: list[float] = []
            for rows, cols in ((li, ri), (ri, li)):
                scores: list[float] = []
                for i in rows:
                    values = [float(pae[i][j]) for j in cols if j < len(pae[i]) and isfinite(float(pae[i][j]))]
                    blocks.extend(values)
                    good = [value for value in values if value < cutoff]
                    if good:
                        scores.append(sum(1.0 / (1.0 + (value / scale) ** 2) for value in good) / len(good))
                if scores:
                    directional.append(max(scores))
            if blocks:
                result[f"{left}:{right}"] = (sum(blocks) / len(blocks), max(directional) if directional else 0.0)
    return result
