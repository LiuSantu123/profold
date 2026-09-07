#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 3 ]]; then
  echo "Usage: $0 INPUT_FASTA RAW_CSV EXTRACTED_CSV" >&2
  exit 2
fi

INPUT_FASTA="$1"
RAW_CSV="$2"
EXTRACTED_CSV="$3"
PYTHON_BIN="${NETSOLP_PYTHON:-/xcfhome/yhliu/002_software/001_conda/envs/screening/bin/python}"
PREDICT_SCRIPT="${NETSOLP_SCRIPT:-/xcfhome/yhliu/002_software/088_NetSolP/PredictionServer/predict_v2.py}"
EXTRACT_SCRIPT="${NETSOLP_EXTRACT_SCRIPT:-/xcfhome/yhliu/003_scripts/01_prediction/extract_netsolp.py}"

[[ -f "$INPUT_FASTA" ]] || { echo "NetSolP FASTA not found: $INPUT_FASTA" >&2; exit 1; }
mkdir -p "$(dirname "$RAW_CSV")" "$(dirname "$EXTRACTED_CSV")"

"$PYTHON_BIN" "$PREDICT_SCRIPT" \
  --FASTA_PATH "$INPUT_FASTA" \
  --OUTPUT_PATH "$RAW_CSV" \
  --MODEL_TYPE "${NETSOLP_MODEL_TYPE:-Both}" \
  --PREDICTION_TYPE "${NETSOLP_PREDICTION_TYPE:-SU}"
"$PYTHON_BIN" "$EXTRACT_SCRIPT" "$RAW_CSV" "$EXTRACTED_CSV"
