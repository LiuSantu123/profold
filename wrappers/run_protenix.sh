#!/usr/bin/env bash
set -euo pipefail
[[ $# -eq 2 ]] || { echo "Usage: $0 INPUT_JSON OUTPUT_DIR" >&2; exit 2; }
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# The v38 Level 2 worker passes a directory so Protenix loads the model once
# for the complete shard.  Keep the historical single-file launcher available
# for older callers, but use the native batch CLI when a directory is given.
if [[ -d "$1" ]]; then
  protenix_bin="${PROTENIX_BIN:-protenix}"
  exec "$protenix_bin" pred \
    --input "$1" --out_dir "$2" \
    --model_name "${PROTENIX_MODEL:-protenix_base_default_v1.0.0}" \
    --seeds "${PROTENIX_SEED:-101}" --cycle "${PROTENIX_CYCLE:-10}" \
    --step "${PROTENIX_STEP:-200}" --sample "${PROTENIX_SAMPLE:-1}" \
    --dtype "${PROTENIX_DTYPE:-bf16}" --use_msa false \
    --need_atom_confidence true
fi

: "${PROTENIX_PYTHON:?Set PROTENIX_PYTHON to the Protenix environment Python}"
: "${PROTENIX_SOURCE:?Set PROTENIX_SOURCE to the matching Protenix source checkout}"
export PYTHONPATH="$PROTENIX_SOURCE${PYTHONPATH:+:$PYTHONPATH}"
exec "$PROTENIX_PYTHON" "${PROTENIX_SCRIPT:-$SCRIPT_DIR/protenix_100.py}" \
  --model_name "${PROTENIX_MODEL:-protenix_base_20250630_v1.0.0}" \
  --seeds "${PROTENIX_SEED:-101}" --dump_dir "$2" --input_json_path "$1" \
  --model.N_cycle 10 --sample_diffusion.N_sample 5 --sample_diffusion.N_step 200 \
  --triangle_attention "${PROTENIX_ATTENTION:-triattention}" \
  --triangle_multiplicative "${PROTENIX_MULTIPLICATIVE:-cuequivariance}" \
  --use_msa false --need_atom_confidence true
