#!/usr/bin/env python3
"""Adapt Level1's multi-record FASTA contract to the local ESMFold2 batch wrapper."""
from __future__ import annotations
import argparse, subprocess, tempfile
from pathlib import Path

def main() -> int:
    ap=argparse.ArgumentParser(); ap.add_argument('-i','--input',type=Path,required=True); ap.add_argument('-o','--output',type=Path,required=True); a=ap.parse_args()
    records=[]; name=None; chunks=[]
    for line in a.input.read_text().splitlines():
        if line.startswith('>'):
            if name is not None: records.append((name,''.join(chunks)))
            name=line[1:].split()[0]; chunks=[]
        elif line.strip(): chunks.append(line.strip())
    if name is not None: records.append((name,''.join(chunks)))
    if not records: raise SystemExit('no FASTA records')
    with tempfile.TemporaryDirectory(prefix='profold_l1_esm_') as tmp:
        inp=Path(tmp)
        for design, seq in records: (inp/f'{design}.fasta').write_text(f'>{design}\n{seq}\n')
        wrapper=Path(__file__).with_name('run_esmfold2_local_batch.sh')
        return subprocess.call(['bash',str(wrapper),'-i',str(inp),'-o',str(a.output)])
if __name__=='__main__': raise SystemExit(main())
