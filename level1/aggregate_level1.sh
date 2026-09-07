#!/usr/bin/env bash
#SBATCH --job-name=v37_l1_agg
#SBATCH --output=level1_aggregate.out
#SBATCH --error=level1_aggregate.err
#SBATCH --cpus-per-task=2
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PYTHON_BIN="${PYTHON_BIN:-python}"
exec "$PYTHON_BIN" "$SCRIPT_DIR/aggregate_level1.py" --run-dir "$1"

