#!/usr/bin/env bash
set -euo pipefail

# Batch adapter for Level 2: convert the shard's FASTA files to a JSON
# directory once, then launch OpenDDE once for the complete directory.  This
# keeps model loading outside the per-design loop.

INPUT=""
OUTPUT=""
while [[ $# -gt 0 ]]; do
  case "$1" in
    -i|--input) INPUT="$2"; shift 2 ;;
    -o|--output|--out_dir) OUTPUT="$2"; shift 2 ;;
    *) echo "unknown option: $1" >&2; exit 2 ;;
  esac
done
[[ -d "$INPUT" && -n "$OUTPUT" ]] || {
  echo "usage: $0 -i FASTA_DIR -o OUTPUT_DIR" >&2
  exit 2
}
mkdir -p "$OUTPUT"
json_dir="$OUTPUT/.opendde_inputs"
mkdir -p "$json_dir"
converter="${OPENDDE_FASTA2JSON:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/fasta_to_opendde_json.py}"
PYTHON_BIN="${PYTHON_BIN:-python}"
for fasta in "$INPUT"/*.fasta; do
  [[ -f "$fasta" ]] || continue
  "$PYTHON_BIN" "$converter" "$fasta" -o "$json_dir/$(basename "${fasta%.fasta}").json"
done

opendde_bin="${OPENDDE_BIN:-opendde}"
# Batch jobs often run without activating a conda environment.  Prefer the
# explicit override, then a command on PATH, and finally the standard local
# OpenDDE environment used by this host.
if [[ -z "${OPENDDE_BIN:-}" ]] && ! command -v "$opendde_bin" >/dev/null 2>&1 \
   && [[ -x /home/amax/.conda/envs/opendde/bin/opendde ]]; then
  opendde_bin=/home/amax/.conda/envs/opendde/bin/opendde
fi
exec "$opendde_bin" pred \
  -i "$json_dir" -o "$OUTPUT" \
  -n "${OPENDDE_MODEL:-opendde_v1}" \
  -s "${OPENDDE_SEEDS:-101}" \
  -c "${OPENDDE_CYCLE:-10}" -p "${OPENDDE_STEP:-200}" \
  -e "${OPENDDE_SAMPLE:-1}" -d "${OPENDDE_DTYPE:-bf16}" \
  --device "${OPENDDE_DEVICE:-cuda}" \
  --use_msa false --need_atom_confidence true
