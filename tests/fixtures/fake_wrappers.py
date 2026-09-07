#!/usr/bin/env python3
from __future__ import annotations

import csv
import json
import sys
from pathlib import Path


def ids_from_fasta(path: Path) -> list[str]:
    return [line[1:].split()[0] for line in path.read_text().splitlines() if line.startswith(">")]


def write_structure(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("data_fake\n")


def main() -> int:
    argv = sys.argv[1:]
    kind = argv[argv.index("--kind") + 1]
    args = argv[: argv.index("--kind")] + argv[argv.index("--kind") + 2 :]

    def value(flag: str) -> Path:
        return Path(args[args.index(flag) + 1])

    if kind == "esmfold":
        input_path, outdir = value("--input"), value("--outdir")
        outdir.mkdir(parents=True, exist_ok=True)
        ids = ids_from_fasta(input_path)
        with (outdir / "esmfold2_summary.csv").open("w", newline="") as handle:
            writer = csv.writer(handle)
            writer.writerow(["design_id", "mean_plddt"])
            for design_id in ids:
                writer.writerow([design_id, "88.0"])
                write_structure(outdir / f"{design_id}.pdb")
        return 0

    if kind == "af3":
        input_dir, outdir = value("--input_dir"), value("--output_dir")
        outdir.mkdir(parents=True, exist_ok=True)
        for source in input_dir.glob("*.json"):
            design_id = source.stem
            (outdir / f"{design_id}_summary_confidences.json").write_text(
                json.dumps({"ranking_score": 0.9, "iptm": 0.8, "ptm": 0.85})
            )
            write_structure(outdir / f"{design_id}.cif")
        return 0

    if kind == "for_wj_cid":
        fasta = value("--fasta")
        json_dir = value("--json-output-dir")
        outdir = value("--af3-output-dir")
        if list(json_dir.glob("*.json")):
            raise SystemExit("unexpected prebuilt AF3 JSON for for_wj_cid")
        json_dir.mkdir(parents=True, exist_ok=True)
        outdir.mkdir(parents=True, exist_ok=True)
        for design_id in ids_from_fasta(fasta):
            for variant, score in (("with_lig", 0.91), ("without_lig", 0.81)):
                prefix = f"{design_id}_af3_{variant}"
                (json_dir / f"{prefix}.json").write_text(json.dumps({"name": prefix}))
                (outdir / f"{prefix}_summary_confidences.json").write_text(
                    json.dumps({"ranking_score": score, "iptm": score - 0.1, "ptm": score - 0.05})
                )
                write_structure(outdir / f"{prefix}.cif")
        return 0

    if kind == "protenix":
        input_path, outdir = Path(args[0]), Path(args[1])
        payload = json.loads(input_path.read_text())
        design_id = payload[0]["name"]
        outdir.mkdir(parents=True, exist_ok=True)
        (outdir / f"{design_id}_summary_confidence_sample_0.json").write_text(
            json.dumps({"plddt": 87.0, "iptm": 0.75})
        )
        write_structure(outdir / f"{design_id}.cif")
        return 0

    if kind == "boltz2":
        input_path, outdir = value("--out_dir"), value("--out_dir")
        input_path = Path(args[args.index("predict") + 1])
        design_id = input_path.stem
        outdir.mkdir(parents=True, exist_ok=True)
        (outdir / f"confidence_{design_id}.json").write_text(
            json.dumps({"plddt": 0.86, "complex_iplddt": 0.74})
        )
        write_structure(outdir / f"{design_id}.cif")
        return 0

    if kind == "opendde":
        input_path, outdir = value("-i"), value("-o")
        design_id = input_path.stem
        outdir.mkdir(parents=True, exist_ok=True)
        (outdir / f"{design_id}_summary_confidences.json").write_text(
            json.dumps({"plddt": 84.0, "iptm": 0.7})
        )
        write_structure(outdir / f"{design_id}.cif")
        return 0

    if kind == "netsolp":
        fasta, raw, extracted = map(Path, args[:3])
        ids = ids_from_fasta(fasta)
        raw.parent.mkdir(parents=True, exist_ok=True)
        with extracted.open("w", newline="") as handle:
            writer = csv.writer(handle)
            writer.writerow(["sid", "predicted_solubility", "predicted_usability"])
            for design_id in ids:
                writer.writerow([design_id, "0.7", "0.8"])
        raw.write_text("fake raw\n")
        return 0

    if kind == "temberture":
        fasta, output = Path(args[0]), value("--output")
        output.parent.mkdir(parents=True, exist_ok=True)
        with output.open("w", newline="") as handle:
            writer = csv.writer(handle)
            writer.writerow(["sid", "tem_TM", "tem_SC"])
            for design_id in ids_from_fasta(fasta):
                writer.writerow([design_id, "62.5", "0.9"])
        return 0

    raise SystemExit(f"unknown kind: {kind}")


if __name__ == "__main__":
    raise SystemExit(main())
