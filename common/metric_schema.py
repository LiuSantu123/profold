#!/usr/bin/env python3
"""Canonical metric names shared by predictors and downstream screens.

Predictor adapters may expose additional model-specific columns, but these
names are stable across AF3, Protenix, Boltz, OpenDDE and ESMFold2.  Per-chain
and per-pair values remain JSON columns so chain IDs are not lost in a flat
table.
"""

from __future__ import annotations

import json
from typing import Mapping


CANONICAL_FIELDS = (
    "design_id", "model", "state", "status",
    "complex_iptm", "complex_ptm", "complex_mean_plddt",
    "complex_mean_pae", "complex_mean_interface_pae", "complex_ipsae",
    "chain_metrics_json", "chain_pair_metrics_json",
    "structure_path", "summary_path", "full_data_path", "full_data_npz_path",
    "full_data_archive_members",
    "reference_structure_path", "usalign_rmsd", "usalign_tm_score",
    "usalign_chain_metrics_json",
    "rosetta_interface_complementarity", "rosetta_sasa",
    "rosetta_delta_sasa", "rosetta_radius_of_gyration", "rosetta_isoelectric_point",
    "geometry_radius_of_gyration", "geometry_sasa", "geometry_delta_sasa",
    "geometry_isoelectric_point", "geometry_source_structure",
    "geometry_isoelectric_point",
    "geometry_interface_pairs_json",
    "geometry_interface_contact_density",
)


def _first(row: Mapping[str, object], *names: str) -> object:
    for name in names:
        value = row.get(name)
        if value not in (None, ""):
            return value
    return ""


def canonical_metric_row(row: Mapping[str, object], model: str | None = None, state: str = "single") -> dict[str, object]:
    """Map a predictor-prefixed row into the stable Level 2 metric schema."""
    prefix = model or str(row.get("model", ""))
    def value(suffix: str, *aliases: str) -> object:
        names = (f"{prefix}_{suffix}",) + tuple(aliases)
        return _first(row, *names)
    output = {field: "" for field in CANONICAL_FIELDS}
    output.update({
        "design_id": row.get("design_id", ""),
        "model": prefix,
        "state": state,
        "status": row.get(f"{prefix}_status", row.get("status", "")),
        "complex_iptm": value("iptm"),
        "complex_ptm": value("ptm"),
        "complex_mean_plddt": value("mean_plddt", f"{prefix}_plddt", f"{prefix}_complex_plddt"),
        "complex_mean_pae": value("mean_pae", f"{prefix}_avg_pae", f"{prefix}_average_pae"),
        "complex_mean_interface_pae": value("mean_interface_pae", f"{prefix}_avg_if_pae"),
        "complex_ipsae": value("ipsae", f"{prefix}_overall_ipsae_mean"),
        "structure_path": value("structure_path"),
        "summary_path": value("summary_json"),
        "full_data_path": value("full_data_json"),
        "full_data_npz_path": value("full_data_npz"),
        "reference_structure_path": row.get("reference_structure_path", ""),
        "usalign_rmsd": row.get(f"{prefix}_usalign_rmsd", ""),
        "usalign_tm_score": row.get(f"{prefix}_usalign_tm_score", ""),
    })
    for field, suffixes in {
        "chain_metrics_json": ("chain_metrics", "chain_ptm", "chain_plddt"),
        "chain_pair_metrics_json": ("chain_pair_metrics", "chain_pair_iptm"),
    }.items():
        raw = value(suffixes[0])
        if raw == "":
            raw = {suffix: value(suffix) for suffix in suffixes if value(suffix) != ""}
        if raw not in ("", {}):
            output[field] = raw if isinstance(raw, str) else json.dumps(raw, ensure_ascii=False, sort_keys=True)
    return output
