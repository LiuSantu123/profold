#!/usr/bin/env bash
#SBATCH --job-name=v38_1_l2
#SBATCH --gres=gpu:1
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.err
set -euo pipefail
# Slurm installations may copy this batch script into a spool directory before
# execution.  Resolve the source tree from an explicit override or the checked
# out ProFold location instead of relying on BASH_SOURCE after that copy.
SCRIPT_DIR="${PROFOLD_ROOT:-/data/codex/tools/profold}/level2"
PYTHON_BIN="${PYTHON_BIN:-python}"
if [[ -n "${PROFOLD_CONDA_ENV:-}" ]] && command -v conda >/dev/null 2>&1; then
  source "$(conda info --base)/etc/profile.d/conda.sh"
  conda activate "$PROFOLD_CONDA_ENV"
  PYTHON_BIN="${PYTHON_BIN:-python}"
fi
exec "$PYTHON_BIN" "$SCRIPT_DIR/run_level2.py" \
  --run-dir "$1" \
  --shard-index "$SLURM_ARRAY_TASK_ID"
