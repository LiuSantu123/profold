#!/usr/bin/env python3
"""Declarative Level 2 predictor registry.

The registry contains only contracts.  Site-specific commands stay in JSON
configuration, while the worker uses these records for input and parser
dispatch.  Adding a predictor therefore does not require editing every
batch loop.
"""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class PredictorSpec:
    name: str
    input_field: str
    input_suffix: str
    parser: str
    batch_placeholder: str


PREDICTORS: dict[str, PredictorSpec] = {
    "protenix": PredictorSpec("protenix", "protenix_input_path", ".json", "summary", "protenix_input_dir"),
    "boltz2": PredictorSpec("boltz2", "boltz_input_path", ".yaml", "summary", "boltz_input_dir"),
    "opendde": PredictorSpec("opendde", "opendde_input_path", ".fasta", "summary", "opendde_input_dir"),
    "esmfold2": PredictorSpec("esmfold2", "esmfold2_input_path", ".fasta", "esmfold2", "esmfold2_input_dir"),
}


def get_predictor(name: str) -> PredictorSpec:
    try:
        return PREDICTORS[name]
    except KeyError as exc:
        raise ValueError(f"unsupported level2 predictor: {name}") from exc
