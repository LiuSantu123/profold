#!/usr/bin/env bash
set -euo pipefail

# The historical wrapper defaults to /pubhome.  Activate the local environment
# first, then delegate so its FASTA/JSON conversion and model flags stay shared.
CONDA_BASE="${CONDA_BASE:-/xcfhome/yhliu/002_software/001_conda}"
OPENDDE_ENV="${OPENDDE_ENV:-opendde}"
HISTORICAL_WRAPPER="${OPENDDE_HISTORICAL_WRAPPER:-/xcfhome/yhliu/002_software/110_opendde/run_opendde.sh}"

[[ -f "$HISTORICAL_WRAPPER" ]] || {
  echo "OpenDDE wrapper not found: $HISTORICAL_WRAPPER" >&2
  exit 1
}

# Ensure the delegated script sees the /xcfhome installation even when the
# submission shell inherited another conda environment.
source "$CONDA_BASE/etc/profile.d/conda.sh"
conda activate "$CONDA_BASE/envs/$OPENDDE_ENV"
export CONDA_BASE OPENDDE_ENV
exec bash "$HISTORICAL_WRAPPER" "$@"
