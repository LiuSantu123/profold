#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PYTHON_BIN="${PYTHON_BIN:-python}"
MANIFEST=""
OUTDIR=""
SHARD_SIZE="32"
CONFIG=""
DEPENDENCY=""
PARTITION=""
TIME_LIMIT=""

usage() {
  echo "Usage: $0 --manifest INPUT.tsv --outdir RUN_DIR [--shard-size N] [--config config.json] [--partition PARTITION] [--time TIME] [--dependency JOBID]"
}
while [[ $# -gt 0 ]]; do
  case "$1" in
    --manifest) MANIFEST="$2"; shift 2 ;;
    --outdir) OUTDIR="$2"; shift 2 ;;
    --shard-size) SHARD_SIZE="$2"; shift 2 ;;
    --config) CONFIG="$2"; shift 2 ;;
    --partition) PARTITION="$2"; shift 2 ;;
    --time) TIME_LIMIT="$2"; shift 2 ;;
    --dependency) DEPENDENCY="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "unknown argument: $1" >&2; usage >&2; exit 2 ;;
  esac
done
[[ -n "$MANIFEST" && -n "$OUTDIR" ]] || { usage >&2; exit 2; }

prepare_args=(--manifest "$MANIFEST" --outdir "$OUTDIR" --shard-size "$SHARD_SIZE")
[[ -n "$CONFIG" ]] && prepare_args+=(--config "$CONFIG")
"$PYTHON_BIN" "$SCRIPT_DIR/prepare_level1.py" "${prepare_args[@]}"
n_shards="$(tr -d '[:space:]' < "$OUTDIR/tasks/n_shards.txt")"
array_args=(--parsable --array="0-$((n_shards - 1))")
[[ -n "$PARTITION" ]] && array_args+=(--partition="$PARTITION")
[[ -n "$TIME_LIMIT" ]] && array_args+=(--time="$TIME_LIMIT")
[[ -n "$DEPENDENCY" ]] && array_args+=(--dependency="afterok:$DEPENDENCY")
array_job="$(sbatch "${array_args[@]}" "$SCRIPT_DIR/run_level1_array.sh" "$OUTDIR")"
aggregate_args=(--parsable --dependency="afterok:$array_job")
[[ -n "$PARTITION" ]] && aggregate_args+=(--partition="$PARTITION")
[[ -n "$TIME_LIMIT" ]] && aggregate_args+=(--time="$TIME_LIMIT")
aggregate_job="$(sbatch "${aggregate_args[@]}" "$SCRIPT_DIR/aggregate_level1.sh" "$OUTDIR")"
echo "level1_array_job=$array_job"
echo "level1_aggregate_job=$aggregate_job"
