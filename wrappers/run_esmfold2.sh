#!/usr/bin/env bash
set -euo pipefail

# Adapter contract for the Level 2 ESMFold2 stage.
#
# ESMFold2 distributions currently expose different launchers.  A site can
# provide an exact command through ESMFOLD2_COMMAND, using {input} and
# {outdir}; otherwise this adapter tries the common `esmfold2 predict` form.
# The model must write at least one PDB containing pLDDT in the B-factor
# column.  The adapter normalizes a single PDB to <design_id>.pdb so the
# ProFold parser can validate it and resume safely.

usage() {
  echo "usage: $0 -i INPUT_FASTA -o OUTPUT_DIR" >&2
  exit 2
}

input=""
outdir=""
while getopts ":i:o:" opt; do
  case "$opt" in
    i) input="$OPTARG" ;;
    o) outdir="$OPTARG" ;;
    *) usage ;;
  esac
done
[[ -n "$input" && -n "$outdir" ]] || usage
[[ -s "$input" ]] || { echo "input FASTA is missing: $input" >&2; exit 2; }
mkdir -p "$outdir"

design_id="$(awk '/^>/{sub(/^>/, ""); split($0, a, /[[:space:]]/); print a[1]; exit}' "$input")"
[[ -n "$design_id" ]] || { echo "could not read design id from FASTA" >&2; exit 2; }

if [[ -n "${ESMFOLD2_COMMAND:-}" ]]; then
  command_string="${ESMFOLD2_COMMAND//\{input\}/$input}"
  command_string="${command_string//\{outdir\}/$outdir}"
  bash -lc "$command_string"
else
  esmfold2_bin="${ESMFOLD2_BIN:-esmfold2}"
  "$esmfold2_bin" predict --input "$input" --output-dir "$outdir"
fi

# Common implementations emit exactly one PDB but use a generic name.  Keep
# all original artifacts and add the stable ProFold filename only when it is
# unambiguous.  Do not overwrite an existing normalized result.
if [[ ! -s "$outdir/$design_id.pdb" ]]; then
  mapfile -t pdbs < <(find "$outdir" -maxdepth 2 -type f -name '*.pdb' -print | sort)
  if [[ "${#pdbs[@]}" -eq 1 ]]; then
    cp -n "${pdbs[0]}" "$outdir/$design_id.pdb"
  fi
fi

[[ -s "$outdir/$design_id.pdb" ]] || {
  echo "ESMFold2 completed without an unambiguous PDB for $design_id" >&2
  exit 3
}
