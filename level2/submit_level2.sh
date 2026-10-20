#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PYTHON_BIN="${PYTHON_BIN:-python}"
MANIFEST=""
OUTDIR=""
SHARD_SIZE="8"
CONFIG=""
DEPENDENCY=""
PARTITION=""
TIME_LIMIT=""
CID_MODE=""
CONDA_ENV=""
DELIVERY_DIR=""

usage() {
  echo "Usage: $0 --manifest INPUT.tsv --outdir RUN_DIR [--shard-size N] [--config config.json] [--conda-env ENV] [--delivery-dir DIR] [--cid] [--partition PARTITION] [--time TIME] [--dependency JOBID]"
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
    --cid) CID_MODE="--cid"; shift ;;
    --conda-env) CONDA_ENV="$2"; shift 2 ;;
    --delivery-dir) DELIVERY_DIR="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "unknown argument: $1" >&2; usage >&2; exit 2 ;;
  esac
done
[[ -n "$MANIFEST" && -n "$OUTDIR" ]] || { usage >&2; exit 2; }
[[ -n "$DELIVERY_DIR" ]] || DELIVERY_DIR="$OUTDIR/review"

prepare_args=(--manifest "$MANIFEST" --outdir "$OUTDIR" --shard-size "$SHARD_SIZE")
[[ -n "$CONFIG" ]] && prepare_args+=(--config "$CONFIG")
[[ -n "$CID_MODE" ]] && prepare_args+=("$CID_MODE")
"$PYTHON_BIN" "$SCRIPT_DIR/prepare_level2.py" "${prepare_args[@]}"
n_shards="$(tr -d '[:space:]' < "$OUTDIR/tasks/n_shards.txt")"
array_args=(--parsable --array="0-$((n_shards - 1))")
[[ -n "$CONDA_ENV" ]] && array_args+=(--export="ALL,PROFOLD_CONDA_ENV=$CONDA_ENV")
[[ -n "$PARTITION" ]] && array_args+=(--partition="$PARTITION")
[[ -n "$TIME_LIMIT" ]] && array_args+=(--time="$TIME_LIMIT")
[[ -n "$DEPENDENCY" ]] && array_args+=(--dependency="afterok:$DEPENDENCY")
array_job="$(sbatch "${array_args[@]}" "$SCRIPT_DIR/run_level2_array.sh" "$OUTDIR")"
aggregate_args=(--parsable --dependency="afterok:$array_job")
[[ -n "$CONDA_ENV" ]] && aggregate_args+=(--export="ALL,PROFOLD_CONDA_ENV=$CONDA_ENV")
[[ -n "$PARTITION" ]] && aggregate_args+=(--partition="$PARTITION")
[[ -n "$TIME_LIMIT" ]] && aggregate_args+=(--time="$TIME_LIMIT")
aggregate_job="$(sbatch "${aggregate_args[@]}" "$SCRIPT_DIR/aggregate_level2.sh" "$OUTDIR")"
closeout_args=(--parsable --dependency="afterok:$aggregate_job")
[[ -n "$PARTITION" ]] && closeout_args+=(--partition="$PARTITION")
[[ -n "$TIME_LIMIT" ]] && closeout_args+=(--time="$TIME_LIMIT")
[[ -n "$CONDA_ENV" ]] && closeout_args+=(--export="ALL,PROFOLD_CONDA_ENV=$CONDA_ENV")
closeout_job="$(sbatch "${closeout_args[@]}" "$SCRIPT_DIR/closeout_level2.sh" "$OUTDIR" "$DELIVERY_DIR")"
echo "level2_array_job=$array_job"
echo "level2_aggregate_job=$aggregate_job"
echo "level2_closeout_job=$closeout_job"
echo "level2_delivery_dir=$DELIVERY_DIR"
