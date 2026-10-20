#!/usr/bin/env bash
set -euo pipefail
input=""; output=""
while [[ $# -gt 0 ]]; do
  case "$1" in
    -i|--input) input="$2"; shift 2;;
    -o|--output|--out_dir) output="$2"; shift 2;;
    *) echo "unknown option: $1" >&2; exit 2;;
  esac
done
[[ -f "$input" && -n "$output" ]] || { echo "usage: $0 -i FASTA -o OUTPUT_DIR" >&2; exit 2; }
exec python "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/esmfold2_level1_adapter.py" -i "$input" -o "$output"
