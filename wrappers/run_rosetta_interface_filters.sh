#!/usr/bin/env bash
set -euo pipefail
ROOT="${PROFOLD_ROOT:-/data/codex/tools/profold}"
ROSETTA_BIN="${ROSETTA_INTERFACE_ANALYZER:-/nvme_af3/softwares/rosetta3.9/rosetta_2016.20/main/source/build/src/release/linux/3.10/64/x86/gcc/4.8/static/InterfaceAnalyzer.static.linuxgccrelease}"
exec python "$ROOT/tools/rosetta/run_interface_filters.py" --input "$2" --output "$4" --executable "$ROSETTA_BIN"
