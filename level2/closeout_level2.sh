#!/usr/bin/env bash
#SBATCH --job-name=v38_1_l2_closeout
#SBATCH --cpus-per-task=1
#SBATCH --output=level2_closeout_%j.out
#SBATCH --error=level2_closeout_%j.err
set -euo pipefail

SCRIPT_DIR="${PROFOLD_ROOT:-/data/codex/tools/profold}/level2"
PYTHON_BIN="${PYTHON_BIN:-python}"
if [[ -n "${PROFOLD_CONDA_ENV:-}" ]] && command -v conda >/dev/null 2>&1; then
  source "$(conda info --base)/etc/profile.d/conda.sh"
  conda activate "$PROFOLD_CONDA_ENV"
fi
RUN_DIR="${1:?run directory is required}"
DELIVERY_DIR="${2:-$RUN_DIR/review}"
exec "$PYTHON_BIN" "$SCRIPT_DIR/compact_level2.py" \
  --run-dir "$RUN_DIR" \
  --out-dir "$DELIVERY_DIR" \
  --delivery-only \
  --archive
