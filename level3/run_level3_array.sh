#!/usr/bin/env bash
#SBATCH --job-name=v37_l3
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=4
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.err
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PYTHON_BIN="${PYTHON_BIN:-python}"
exec "$PYTHON_BIN" "$SCRIPT_DIR/run_level3.py" \
  --run-dir "$1" \
  --shard-index "$SLURM_ARRAY_TASK_ID"
