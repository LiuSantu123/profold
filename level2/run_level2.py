#!/usr/bin/env python3
"""Run one Level 2 shard with three independently resumable model stages."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

V37 = Path(__file__).resolve().parents[1]
if str(V37) not in sys.path:
    sys.path.insert(0, str(V37))

from common.input_builder import (  # noqa: E402
    build_boltz_yaml,
    build_opendde_fasta,
    build_protenix_json,
    status_from_models,
)
from common.model_parsers import parse_summary_model  # noqa: E402
from common.pipeline import load_config, load_shard, write_part  # noqa: E402
from common.runner import command_template, configured_command, run_logged  # noqa: E402


PROTENIX_WRAPPER = Path("/xcfhome/yhliu/003_scripts/01_prediction/protenix_100_nomsa.sh")
OPENDDE_WRAPPER = V37 / "wrappers" / "run_opendde.sh"
LEVEL2_MODELS = ("protenix", "boltz2", "opendde")


def _configured(config: dict[str, object], name: str, context: dict[str, object]) -> list[str] | None:
    command = configured_command(config, name, context)
    if command is None and name == "boltz2":
        command = configured_command(config, "boltz", context)
    return command


def _default_command(
    config: dict[str, object], name: str, context: dict[str, object]
) -> list[str] | None:
    configured = _configured(config, name, context)
    if configured is not None:
        return configured
    if name == "protenix":
        wrapper = str(config.get("protenix_wrapper", PROTENIX_WRAPPER))
        return command_template(
            ["bash", wrapper, "{protenix_input}", "{protenix_outdir}"], context
        )
    if name == "opendde":
        wrapper = str(config.get("opendde_wrapper", OPENDDE_WRAPPER))
        return command_template(
            ["bash", wrapper, "-i", "{opendde_input}", "-o", "{opendde_outdir}"], context
        )
    if name == "boltz2":
        boltz_bin = str(config.get("boltz_bin", "boltz"))
        command = [
            boltz_bin,
            "predict",
            "{boltz_input}",
            "--out_dir",
            "{boltz_outdir}",
            "--write_full_pae",
            "--write_full_pde",
        ]
        if config.get("boltz_no_kernels", False):
            command.append("--no_kernels")
        return command_template(command, context)
    raise ValueError(f"unsupported level2 model: {name}")


def _build_input(row: dict[str, str], model: str, input_dir: Path) -> str:
    design_id = row["design_id"]
    if model == "protenix":
        return build_protenix_json(row, input_dir / f"{design_id}.json")
    if model == "boltz2":
        return build_boltz_yaml(row, input_dir / f"{design_id}.yaml")
    if model == "opendde":
        return build_opendde_fasta(row, input_dir / f"{design_id}.fasta")
    raise ValueError(model)


def _parse_model(model_dir: Path, design_id: str, model: str) -> dict[str, object]:
    return parse_summary_model(model_dir, design_id, model)


def _stage_empty(model: str, status: str, error: str) -> dict[str, object]:
    return {
        f"{model}_status": status,
        f"{model}_error": error,
        f"{model}_log_path": "",
    }


def _run_one(
    row: dict[str, str],
    model: str,
    run_dir: Path,
    shard_dir: Path,
    config: dict[str, object],
    *,
    resume: bool,
) -> dict[str, object]:
    design_id = row["design_id"]
    model_dir = shard_dir / "models" / model / design_id
    input_dir = shard_dir / "inputs" / model
    model_dir.mkdir(parents=True, exist_ok=True)
    input_dir.mkdir(parents=True, exist_ok=True)
    try:
        input_path = _build_input(row, model, input_dir)
    except Exception as exc:
        return _stage_empty(model, "missing", f"input: {exc}")

    parsed = _parse_model(model_dir, design_id, model)
    if resume and parsed.get(f"{model}_status") == "success":
        parsed[f"{model}_reused"] = "1"
        parsed[f"{model}_log_path"] = ""
        parsed[f"{model}_input_path"] = input_path
        return parsed

    log_path = run_dir / "logs" / f"level2_{model}_{design_id}.log"
    context: dict[str, object] = {
        "run_dir": run_dir,
        "shard_dir": shard_dir,
        "design_id": design_id,
        "model": model,
        "input": input_path,
        "outdir": model_dir,
        "model_input": input_path,
        "model_outdir": model_dir,
        "protenix_input": input_path if model == "protenix" else "",
        "protenix_outdir": model_dir if model == "protenix" else "",
        "boltz_input": input_path if model == "boltz2" else "",
        "boltz_outdir": model_dir if model == "boltz2" else "",
        "opendde_input": input_path if model == "opendde" else "",
        "opendde_outdir": model_dir if model == "opendde" else "",
    }
    command = _default_command(config, model, context)
    if command is None:
        return _stage_empty(model, "missing", f"{model} command is not configured")
    rc, tail = run_logged(command, log_path, cwd=shard_dir)
    parsed = _parse_model(model_dir, design_id, model)
    parsed[f"{model}_log_path"] = str(log_path)
    parsed[f"{model}_input_path"] = input_path
    if rc != 0:
        message = f"command returncode={rc}"
        if tail:
            message = f"{message}; {tail}"
        if parsed.get(f"{model}_status") == "success":
            parsed[f"{model}_status"] = "partial"
        elif parsed.get(f"{model}_status") == "missing":
            parsed[f"{model}_status"] = "failed"
        parsed[f"{model}_error"] = "; ".join(
            item for item in (str(parsed.get(f"{model}_error", "")), message) if item
        )
    return parsed


def run_shard(run_dir: str | Path, shard_index: int, *, resume: bool = True) -> list[dict[str, object]]:
    run = Path(run_dir).resolve()
    config = load_config(run)
    rows = load_shard(run, shard_index)
    shard_dir = run / "artifacts" / f"shard_{shard_index:05d}"
    shard_dir.mkdir(parents=True, exist_ok=True)
    models = LEVEL2_MODELS
    output: list[dict[str, object]] = []
    for row in rows:
        result: dict[str, object] = {"design_id": row["design_id"]}
        for model in models:
            result.update(_run_one(row, model, run, shard_dir, config, resume=resume))
        statuses = {model: str(result.get(f"{model}_status", "missing")) for model in models}
        result["status"] = status_from_models(statuses)
        result["error"] = "; ".join(
            str(result.get(f"{model}_error", "")) for model in models if result.get(f"{model}_error", "")
        )
        result["task_index"] = row.get("task_index", "")
        result["shard_index"] = str(shard_index)
        output.append(result)
    write_part(run, shard_index, output)
    return output


def main() -> int:
    parser = argparse.ArgumentParser(description="Run one v37 Level 2 shard")
    parser.add_argument("--run-dir", required=True)
    parser.add_argument("--shard-index", type=int, required=True)
    parser.add_argument("--no-resume", action="store_true")
    args = parser.parse_args()
    run_shard(args.run_dir, args.shard_index, resume=not args.no_resume)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
