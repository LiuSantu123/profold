#!/usr/bin/env bash
#SBATCH --job-name=v37_l3_agg
#SBATCH --cpus-per-task=1
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PYTHON_BIN="${PYTHON_BIN:-python}"
exec "$PYTHON_BIN" "$SCRIPT_DIR/aggregate_level3.py" --run-dir "$1"
