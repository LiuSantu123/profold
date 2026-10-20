#!/usr/bin/env bash
set -euo pipefail

input=""; output=""
while [[ $# -gt 0 ]]; do
  case "$1" in
    -i|--input) input="$2"; shift 2 ;;
    -o|--output|--out_dir) output="$2"; shift 2 ;;
    *) echo "unknown option: $1" >&2; exit 2 ;;
  esac
done
[[ -d "$input" && -n "$output" ]] || { echo "usage: $0 -i INPUT_DIR -o OUTPUT_DIR" >&2; exit 2; }
export PYTHONPATH="${ESM_VENDOR_ROOT:-/nvme_af3/vendor/esm}${PYTHONPATH:+:$PYTHONPATH}"
export HF_HUB_OFFLINE="${HF_HUB_OFFLINE:-1}"
exec "${ESMFOLD2_PYTHON:-/nvme_af3/conda/envs/esmfold2/bin/python}" \
  "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/esmfold2_local_batch.py" \
  -i "$input" -o "$output"
