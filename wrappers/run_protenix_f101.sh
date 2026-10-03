#!/usr/bin/env bash
set -euo pipefail

[[ $# -eq 2 ]] || { echo "Usage: $0 INPUT_JSON OUTPUT_DIR" >&2; exit 2; }
INPUT_JSON="$1"
OUTPUT_DIR="$2"
PROTENIX_PYTHON="${PROTENIX_PYTHON:-/xcfhome/yhliu/002_software/001_conda/envs/protenix/bin/python}"
PROTENIX_SCRIPT="${PROTENIX_SCRIPT:-/xcfhome/yhliu/003_scripts/01_prediction/protenix_100.py}"

[[ -f "$INPUT_JSON" ]] || { echo "Protenix input not found: $INPUT_JSON" >&2; exit 1; }
mkdir -p "$OUTPUT_DIR"
exec "$PROTENIX_PYTHON" "$PROTENIX_SCRIPT" \
  --model_name protenix_base_20250630_v1.0.0 \
  --seeds 101 \
  --dump_dir "$OUTPUT_DIR" \
  --input_json_path "$INPUT_JSON" \
  --model.N_cycle 10 \
  --sample_diffusion.N_sample 5 \
  --sample_diffusion.N_step 200 \
  --triangle_attention triattention \
  --triangle_multiplicative cuequivariance \
  --use_msa false \
  --need_atom_confidence true
