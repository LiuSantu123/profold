#!/usr/bin/env bash
#SBATCH --job-name=v38_1_l2_agg
# USalign is CPU-bound and runs one process per worker.
# Keep this aligned with the default usalign_workers value.
#SBATCH --cpus-per-task=4
set -euo pipefail
SCRIPT_DIR="${PROFOLD_ROOT:-/data/codex/tools/profold}/level2"
PYTHON_BIN="${PYTHON_BIN:-python}"
if [[ -n "${PROFOLD_CONDA_ENV:-}" ]] && command -v conda >/dev/null 2>&1; then
  source "$(conda info --base)/etc/profile.d/conda.sh"
  conda activate "$PROFOLD_CONDA_ENV"
fi
exec "$PYTHON_BIN" "$SCRIPT_DIR/aggregate_level2.py" --run-dir "$1"
