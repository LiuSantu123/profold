#!/usr/bin/env python3
"""Convert FASTA to OpenDDE input JSON.

支持的多链格式 (参考 examples/dimer.fasta):
  >jobname
  SEQUENCE_CHAIN_A:SEQUENCE_CHAIN_B

- 一条 >record = 一个 job; record 内用 ':' 分隔多条链
- 多条 >record = 多个 job (批量预测)
- 链类型默认 protein; 纯 ACGT 的链可选识别为 dna (加 --auto-type)
- 配体通过 --ligand CCD_CODE 添加 (自动加 CCD_ 前缀), 可多次指定

用法:
  python fasta_to_opendde_json.py input.fasta -o out.json
  python fasta_to_opendde_json.py dimer.fasta -o dimer.json
  python fasta_to_opendde_json.py prot.fasta -o out.json --ligand ATP --ligand MG
  python fasta_to_opendde_json.py prot.fasta -o out.json --ligand 6OI --ligand-count 2
  python fasta_to_opendde_json.py prot.fasta -o out.json --smiles "CC(=O)C"
"""
import argparse
import json
import sys


def parse_fasta(path):
    jobs = []
    name = None
    seq_lines = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    jobs.append((name, "".join(seq_lines)))
                header = line[1:].strip()
                # job name 取第一个 token, 去掉 | 后注释
                name = (header.split("|")[0].split()[0] or "job")
                seq_lines = []
            else:
                seq_lines.append(line)
        if name is not None:
            jobs.append((name, "".join(seq_lines)))
    return jobs


def guess_type(seq, auto_type):
    if not auto_type:
        return "protein"
    s = set(seq.upper())
    if s <= set("ACGT"):
        return "dna"
    if s <= set("ACGU"):
        return "rna"
    return "protein"


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("fasta", help="Input FASTA file")
    ap.add_argument("-o", "--output", required=True, help="Output JSON path")
    ap.add_argument("--ligand", action="append", default=[],
                    help="Ligand CCD code (e.g. ATP, MG, 6OI); auto-prefixed with CCD_. Repeatable.")
    ap.add_argument("--ligand-count", type=int, default=1, help="Count for each --ligand (default 1)")
    ap.add_argument("--smiles", action="append", default=[],
                    help="Ligand as SMILES string. Repeatable.")
    ap.add_argument("--smiles-count", type=int, default=1, help="Count for each --smiles (default 1)")
    ap.add_argument("--auto-type", action="store_true",
                    help="Auto-detect DNA/RNA chains (pure ACGT -> dna, ACGU -> rna). Default: all protein.")
    args = ap.parse_args()

    jobs = parse_fasta(args.fasta)
    if not jobs:
        sys.exit("Error: no sequences found in FASTA")

    out = []
    for name, seq in jobs:
        chains = [c.strip() for c in seq.split(":") if c.strip()]
        if not chains:
            print(f"[WARN] job '{name}' has no sequence, skipped", file=sys.stderr)
            continue
        sequences = []
        cid = ord("A")
        for ch in chains:
            ctype = guess_type(ch, args.auto_type)
            if ctype == "protein":
                sequences.append({"proteinChain": {"sequence": ch, "count": 1, "id": [chr(cid)]}})
            elif ctype == "dna":
                sequences.append({"dnaSequence": {"sequence": ch, "count": 1, "id": [chr(cid)]}})
            elif ctype == "rna":
                sequences.append({"rnaSequence": {"sequence": ch, "count": 1, "id": [chr(cid)]}})
            cid += 1
        for lig in args.ligand:
            code = lig if lig.startswith("CCD_") else "CCD_" + lig
            sequences.append({"ligand": {"ligand": code, "count": args.ligand_count}})
        for smi in args.smiles:
            sequences.append({"ligand": {"smiles": smi, "count": args.smiles_count}})
        out.append({"name": name, "covalent_bonds": [], "sequences": sequences})

    with open(args.output, "w") as f:
        json.dump(out, f, indent=2)
    print(f"Wrote {len(out)} job(s) -> {args.output}")
    for j in out:
        n_prot = sum(1 for s in j["sequences"] if "proteinChain" in s)
        n_dna = sum(1 for s in j["sequences"] if "dnaSequence" in s)
        n_rna = sum(1 for s in j["sequences"] if "rnaSequence" in s)
        n_lig = sum(1 for s in j["sequences"] if "ligand" in s)
        print(f"  - {j['name']}: protein={n_prot} dna={n_dna} rna={n_rna} ligand={n_lig}")


if __name__ == "__main__":
    main()
