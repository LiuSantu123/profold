#!/usr/bin/env bash
#SBATCH --job-name=v37_l1
#SBATCH --output=level1_%A_%a.out
#SBATCH --error=level1_%A_%a.err
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PYTHON_BIN="${PYTHON_BIN:-python}"
: "${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID is required}"
exec "$PYTHON_BIN" "$SCRIPT_DIR/run_level1.py" \
  --run-dir "$1" \
  --shard-index "$SLURM_ARRAY_TASK_ID"

