#!/usr/bin/env python3
"""Run ESMFold2 on a directory of ProFold FASTA inputs.

The released ESMFold2 checkpoint keeps the ESMC-6B backbone in a separate
Hugging Face repository.  Loading by repository name can trigger a network
download even when the six local shards already exist, so this adapter loads
the trunk and the local ESMC snapshot explicitly in offline mode.
"""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path


def read_fasta(path: Path) -> tuple[str, list[tuple[str, str]]]:
    lines = path.read_text(encoding="utf-8").splitlines()
    headers = [line[1:].split()[0] for line in lines if line.startswith(">")]
    if not headers:
        raise ValueError(f"FASTA has no header: {path}")
    name = headers[0]
    sequence = "".join(line.strip() for line in lines if line and not line.startswith(">"))
    chains = [(chr(65 + index), value) for index, value in enumerate(sequence.replace("|", ":").split(":")) if value]
    if not chains:
        raise ValueError(f"FASTA has no sequence: {path}")
    return name, chains


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input-dir", required=True, type=Path)
    parser.add_argument("-o", "--output-dir", required=True, type=Path)
    args = parser.parse_args()
    model_dir = Path(os.environ.get("ESMFOLD2_MODEL_DIR", "/nvme_af3/model/checkpoint/esmfold2/ESMFold2"))
    esmc_dir = Path(os.environ.get(
        "ESMC_MODEL_DIR",
        "/nvme_af3/model/checkpoint/esmfold2/hf_cache/hub/models--biohub--ESMC-6B/snapshots/45b0fa5d7fb06faefbd5e3b89bdcef35d564e79a",
    ))
    if not (model_dir / "config.json").is_file():
        raise SystemExit(f"ESMFold2 model config missing: {model_dir}")
    if not (esmc_dir / "model.safetensors.index.json").is_file():
        raise SystemExit(f"local ESMC snapshot missing: {esmc_dir}")
    inputs = sorted(args.input_dir.glob("*.fasta"))
    if not inputs:
        raise SystemExit(f"no FASTA inputs under {args.input_dir}")

    os.environ.setdefault("HF_HUB_OFFLINE", "1")
    os.environ.setdefault("HF_HOME", "/nvme_af3/model/checkpoint/esmfold2/hf_cache")
    from esm.models.esmfold2 import ESMFold2InputBuilder, EsmFold2Model, ProteinInput, StructurePredictionInput
    import torch

    if not torch.cuda.is_available():
        raise SystemExit("ESMFold2 requires CUDA")
    model = EsmFold2Model.from_pretrained(model_dir, device="cuda", load_esmc=False).eval()
    model.load_esmc(esmc_dir, precision="bf16")
    args.output_dir.mkdir(parents=True, exist_ok=True)
    for fasta in inputs:
        # The manifest filename is the stable design identity.  FASTA headers
        # may come from a shared template (especially for apo/holo inputs), so
        # do not let a template header change the output filename.
        _header_name, chains = read_fasta(fasta)
        design_id = fasta.stem
        structure_input = StructurePredictionInput(
            sequences=[ProteinInput(id=chain_id, sequence=sequence) for chain_id, sequence in chains]
        )
        with torch.inference_mode():
            result = ESMFold2InputBuilder().fold(
                model,
                structure_input,
                num_loops=int(os.environ.get("ESMFOLD2_NUM_LOOPS", "20")),
                num_sampling_steps=int(os.environ.get("ESMFOLD2_NUM_SAMPLING_STEPS", "100")),
                num_diffusion_samples=1,
                seed=0,
            )
        (args.output_dir / f"{design_id}.pdb").write_text(result.complex.to_protein_complex().to_pdb_string())
        (args.output_dir / f"{design_id}.cif").write_text(result.complex.to_mmcif())
        metrics = {
            "design_id": design_id,
            "plddt": 100.0 * float(result.plddt.mean()),
            "ptm": float(result.ptm),
            "sequence_length": sum(len(sequence) for _, sequence in chains),
            "model": "ESMFold2",
            "local": True,
        }
        (args.output_dir / f"{design_id}_confidence.json").write_text(json.dumps(metrics, indent=2) + "\n")
        print(json.dumps(metrics), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
