#!/usr/bin/env bash
# OpenDDE 推理包装脚本
# 自动激活 conda 环境、设置 OPENDDE_ROOT_DIR、调用 opendde pred
# 输入支持 .json (直接用) 或 .fasta/.fa (自动转换为 JSON)
#
# 用法:
#   ./run_opendde.sh -i examples/tiny.json -o ./out
#   ./run_opendde.sh -i prot.fasta -o ./out --ligand ATP
#   ./run_opendde.sh -i dimer.fasta -o ./out -s 101,202 -e 5
#   ./run_opendde.sh -i input.json -o ./out --msa --template   # 完整模式(需 search_database)
set -euo pipefail

# ===== 路径配置 (可被环境变量覆盖) =====
CONDA_BASE="${CONDA_BASE:-/pubhome/yhliu/002_software/001_conda}"
OPENDDE_ENV="${OPENDDE_ENV:-opendde}"
OPENDDE_ROOT="${OPENDDE_ROOT_DIR:-/xcfhome/yhliu/002_software/110_opendde/opendde_data}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
FASTA2JSON="${SCRIPT_DIR}/fasta_to_opendde_json.py"

# ===== 推理默认参数 =====
INPUT=""
OUT_DIR="./opendde_output"
SEEDS="101"
CYCLE=10
STEP=200
SAMPLE=1
DTYPE="bf16"
MODEL="opendde_v1"
USE_MSA=false
USE_TEMPLATE=false
USE_RNA_MSA=false
LIGANDS=()
LIGAND_COUNT=1

usage() {
  cat <<EOF
Usage: $(basename "$0") -i INPUT [-o OUT_DIR] [options]

Required:
  -i, --input FILE       Input JSON or FASTA (.json/.fasta/.fa)

Options:
  -o, --out_dir DIR      Output directory (default: ./opendde_output)
  -s, --seeds STR        Comma-separated seeds (default: 101)
  -c, --cycle N          Pairformer cycles (default: 10)
  -p, --step N           Diffusion steps (default: 200)
  -e, --sample N         Number of samples (default: 1)
  -d, --dtype STR        fp32 | bf16 (default: bf16)
  -n, --model_name STR   Model checkpoint name (default: opendde_v1)
  --msa                  Enable MSA (needs search_database; default off)
  --template             Enable templates (needs search_database; default off)
  --rna_msa              Enable RNA MSA (needs search_database; default off)
  --ligand CCD           Add ligand by CCD code (FASTA mode only; repeatable)
  --ligand-count N       Count for each --ligand (default: 1)
  -h, --help             Show this help

Default = fast no-MSA mode (suitable for monomer / quick screening).
EOF
}

while [ $# -gt 0 ]; do
  case "$1" in
    -i|--input) INPUT="$2"; shift 2;;
    -o|--out_dir) OUT_DIR="$2"; shift 2;;
    -s|--seeds) SEEDS="$2"; shift 2;;
    -c|--cycle) CYCLE="$2"; shift 2;;
    -p|--step) STEP="$2"; shift 2;;
    -e|--sample) SAMPLE="$2"; shift 2;;
    -d|--dtype) DTYPE="$2"; shift 2;;
    -n|--model_name) MODEL="$2"; shift 2;;
    --msa) USE_MSA=true; shift;;
    --template) USE_TEMPLATE=true; shift;;
    --rna_msa) USE_RNA_MSA=true; shift;;
    --ligand) LIGANDS+=("$2"); shift 2;;
    --ligand-count) LIGAND_COUNT="$2"; shift 2;;
    -h|--help) usage; exit 0;;
    *) echo "Unknown option: $1" >&2; usage; exit 1;;
  esac
done

[ -z "$INPUT" ] && { echo "Error: -i INPUT is required" >&2; usage; exit 1; }
[ -f "$INPUT" ] || { echo "Error: input not found: $INPUT" >&2; exit 1; }

# ===== 激活 conda 环境 =====
if [ -z "${CONDA_PREFIX:-}" ] || [ "$(basename "$CONDA_PREFIX")" != "$OPENDDE_ENV" ]; then
  # shellcheck disable=SC1091
  source "${CONDA_BASE}/etc/profile.d/conda.sh"
  conda activate "$OPENDDE_ENV"
fi
export OPENDDE_ROOT_DIR="$OPENDDE_ROOT"

# ===== FASTA -> JSON (如需要) =====
ext="${INPUT##*.}"
case "$ext" in
  fasta|fa|FASTA|FA)
    JSON_INPUT="${OUT_DIR%/}/_input.json"
    mkdir -p "${OUT_DIR%/}"
    lig_args=()
    for l in "${LIGANDS[@]:-}"; do [ -n "$l" ] && lig_args+=(--ligand "$l"); done
    [ "${#lig_args[@]}" -gt 0 ] && lig_args+=(--ligand-count "$LIGAND_COUNT")
    echo "[INFO] Converting FASTA -> JSON: $INPUT -> $JSON_INPUT"
    python "$FASTA2JSON" "$INPUT" -o "$JSON_INPUT" "${lig_args[@]}"
    INPUT="$JSON_INPUT"
    ;;
  json|JSON) ;;
  *) echo "Error: unsupported input extension: $ext (use .json/.fasta/.fa)" >&2; exit 1;;
esac

# ===== 推理 =====
echo "=============================================="
echo " OpenDDE prediction"
echo "  input    : $INPUT"
echo "  out_dir  : $OUT_DIR"
echo "  model    : $MODEL  dtype: $DTYPE"
echo "  seeds    : $SEEDS  sample: $SAMPLE"
echo "  cycle    : $CYCLE  step: $STEP"
echo "  msa      : $USE_MSA  template: $USE_TEMPLATE  rna_msa: $USE_RNA_MSA"
echo "  root_dir : $OPENDDE_ROOT_DIR"
echo "  start    : $(date)"
echo "=============================================="

opendde pred \
  -i "$INPUT" \
  -o "$OUT_DIR" \
  -n "$MODEL" \
  -s "$SEEDS" \
  -c "$CYCLE" \
  -p "$STEP" \
  -e "$SAMPLE" \
  -d "$DTYPE" \
  --use_msa "$USE_MSA" \
  --use_template "$USE_TEMPLATE" \
  --use_rna_msa "$USE_RNA_MSA"

echo "=============================================="
echo " Done: $(date)"
echo " Output: $OUT_DIR"
find "$OUT_DIR" -name "*.cif" -o -name "*confidence*.json" 2>/dev/null | sort
echo "=============================================="
