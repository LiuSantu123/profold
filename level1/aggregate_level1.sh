#!/usr/bin/env bash
#SBATCH --job-name=v37_l1_agg
#SBATCH --output=level1_aggregate.out
#SBATCH --error=level1_aggregate.err
#SBATCH --cpus-per-task=2
set -euo pipefail

SCRIPT_DIR="${PROFOLD_ROOT:-/data/codex/tools/profold}/level1"
PYTHON_BIN="${PYTHON_BIN:-python}"
if [[ -n "${PROFOLD_CONDA_ENV:-}" ]] && command -v conda >/dev/null 2>&1; then
  source "$(conda info --base)/etc/profile.d/conda.sh"
  conda activate "$PROFOLD_CONDA_ENV"
fi
exec "$PYTHON_BIN" "$SCRIPT_DIR/aggregate_level1.py" --run-dir "$1"
