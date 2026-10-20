#!/usr/bin/env python3
"""Public CSV-to-Slurm entry point for ProFold campaigns.

This keeps user-facing input independent from the internal TSV and shard
layout.  ``--dry-run`` prepares and prints the exact submission command, which
is useful on login nodes or while GPUs are occupied.
"""

from __future__ import annotations

import argparse
import json
import shlex
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from common.sequence_csv import write_manifest_from_csv


def _normalize_config_paths(config: dict[str, object], base: Path) -> None:
    """Make site paths independent of the directory from which sbatch runs."""

    keys = (
        "protenix_wrapper", "opendde_wrapper", "esmfold2_wrapper", "boltz_bin",
        "esmfold2_bin", "usalign_path", "rosetta_filter_script",
    )
    for key in keys:
        value = config.get(key)
        if isinstance(value, str) and value and not Path(value).expanduser().is_absolute():
            config[key] = str((base / value).resolve())
    commands = config.get("commands")
    if isinstance(commands, dict):
        for name, raw in list(commands.items()):
            tokens = list(raw) if isinstance(raw, (list, tuple)) else None
            if tokens is None:
                continue
            normalized: list[object] = []
            for token in tokens:
                value = str(token)
                candidate = (base / value).resolve() if value.startswith("./") else None
                normalized.append(str(candidate) if candidate and candidate.is_file() else token)
            commands[name] = normalized


def main() -> int:
    parser = argparse.ArgumentParser(description="Submit a ProFold campaign from sequence CSV")
    parser.add_argument("--sequence-csv", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--config", required=True, help="JSON with model commands and conda/script paths")
    parser.add_argument("--conda-env", help="override the controller/worker conda environment")
    parser.add_argument("--run-script", help="override the level submission script")
    parser.add_argument("--level", choices=("level1", "level2"), default="level2")
    parser.add_argument("--shard-size", type=int, default=8)
    parser.add_argument("--partition")
    parser.add_argument("--time")
    parser.add_argument("--delivery-dir")
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()
    config_path = Path(args.config).expanduser().resolve()
    config = json.loads(config_path.read_text(encoding="utf-8"))
    if not isinstance(config, dict):
        raise SystemExit("--config must contain a JSON object")
    _normalize_config_paths(config, config_path.parent)
    out = Path(args.outdir).expanduser().resolve()
    if args.conda_env:
        config["conda_env"] = args.conda_env
    if args.run_script:
        config["run_script"] = str(Path(args.run_script).expanduser().resolve())
    manifest = out / "input_manifest.tsv"
    write_manifest_from_csv(args.sequence_csv, manifest)
    prepare = ROOT / args.level / f"prepare_{args.level}.py"
    submit = Path(str(config.get("run_script", ROOT / args.level / f"submit_{args.level}.sh"))).expanduser()
    if not submit.is_absolute():
        submit = (config_path.parent / submit).resolve()
    else:
        submit = submit.resolve()
    if not submit.is_file():
        raise SystemExit(f"run script not found: {submit}")
    # Keep config metadata self-contained in the run directory.  The Slurm
    # wrapper consumes this file and remains compatible with direct TSV users.
    run_config = out / "campaign_config.json"
    run_config.parent.mkdir(parents=True, exist_ok=True)
    run_config.write_text(json.dumps(config, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    command = ["bash", str(submit), "--manifest", str(manifest), "--outdir", str(out / "run"),
               "--shard-size", str(args.shard_size), "--config", str(run_config)]
    if config.get("conda_env"):
        command += ["--conda-env", str(config["conda_env"])]
    if args.partition:
        command += ["--partition", args.partition]
    if args.time:
        command += ["--time", args.time]
    if args.delivery_dir:
        command += ["--delivery-dir", str(Path(args.delivery_dir).expanduser().resolve())]
    print("prepared_manifest=" + str(manifest))
    print("submit_command=" + shlex.join(command))
    if args.dry_run:
        return 0
    completed = subprocess.run(command, check=False)
    return completed.returncode


if __name__ == "__main__":
    raise SystemExit(main())
