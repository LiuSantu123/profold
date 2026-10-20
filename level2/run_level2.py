#!/usr/bin/env python3
"""Run one Level 2 shard with independently resumable model stages."""

from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

V37 = Path(__file__).resolve().parents[1]
if str(V37) not in sys.path:
    sys.path.insert(0, str(V37))

from common.input_builder import (  # noqa: E402
    build_boltz_yaml,
    build_esmfold2_fasta,
    build_opendde_fasta,
    build_protenix_json,
    status_from_models,
)
from common.model_parsers import parse_esmfold, parse_summary_model  # noqa: E402
from common.predictor_registry import PREDICTORS, get_predictor  # noqa: E402
from common.pipeline import load_config, load_shard, write_part  # noqa: E402
from common.runner import atomic_json, command_template, configured_command, run_logged  # noqa: E402
from tools.monitors.esmfold2_monitor import scan as scan_esmfold2, write_csv as write_esmfold2_csv  # noqa: E402


PROTENIX_WRAPPER = Path("/xcfhome/yhliu/003_scripts/01_prediction/protenix_100_nomsa.sh")
OPENDDE_WRAPPER = V37 / "wrappers" / "run_opendde.sh"
ESMFOLD2_WRAPPER = V37 / "wrappers" / "run_esmfold2.sh"
LEVEL2_MODELS = tuple(PREDICTORS)
CID_STATES = ("apo", "holo")


def configured_models(config: dict[str, object]) -> tuple[str, ...]:
    """Return the enabled Level 2 predictors, preserving the stable default."""

    raw = config.get("level2_models")
    if raw in (None, ""):
        return LEVEL2_MODELS
    if isinstance(raw, str):
        names = tuple(item.strip() for item in raw.split(",") if item.strip())
    elif isinstance(raw, (list, tuple)):
        names = tuple(str(item).strip() for item in raw if str(item).strip())
    else:
        raise ValueError("level2_models must be a comma-separated string or list")
    unknown = sorted(set(names) - set(LEVEL2_MODELS))
    if unknown or not names:
        raise ValueError(f"unsupported or empty level2_models: {', '.join(unknown) or 'empty'}")
    return names


def _prefix(model: str, state: str) -> str:
    """Keep v38 single-state columns compatible; CID adds state namespace."""
    return model if state == "single" else f"{model}_{state}"


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
            ["bash", wrapper, "{protenix_input_dir}", "{protenix_outdir}"], context
        )
    if name == "opendde":
        wrapper = str(config.get("opendde_wrapper", OPENDDE_WRAPPER))
        return command_template(
            ["bash", wrapper, "-i", "{opendde_input_dir}", "-o", "{opendde_outdir}"], context
        )
    if name == "boltz2":
        boltz_bin = str(config.get("boltz_bin", "boltz"))
        command = [
            boltz_bin,
            "predict",
            "{boltz_input_dir}",
            "--out_dir",
            "{boltz_outdir}",
            "--write_full_pae",
            "--write_full_pde",
        ]
        if config.get("boltz_no_kernels", False):
            command.append("--no_kernels")
        return command_template(command, context)
    if name == "esmfold2":
        wrapper = str(config.get("esmfold2_wrapper", ESMFOLD2_WRAPPER))
        command = ["bash", wrapper, "-i", "{esmfold2_input_dir}", "-o", "{esmfold2_outdir}"]
        if config.get("esmfold2_bin"):
            command = ["env", f"ESMFOLD2_BIN={config['esmfold2_bin']}"] + command
        return command_template(command, context)
    raise ValueError(f"unsupported level2 model: {name}")


def _build_input(row: dict[str, str], model: str, input_dir: Path, state: str) -> str:
    design_id = row["design_id"]
    spec = get_predictor(model)
    normal_field = spec.input_field
    # CID manifests may provide apo_*/holo_* inputs.  If a state-specific
    # input is absent, fall back to the normal input or generated sequence;
    # this keeps single-state manifests valid while making explicit ligand
    # inputs available for real apo/holo screening.
    state_field = f"{state}_{normal_field}"
    input_field = state_field if row.get(state_field) else normal_field
    if spec.parser == "summary" and model == "protenix":
        return build_protenix_json(row, input_dir / f"{design_id}.json", input_field=input_field)
    if spec.parser == "summary" and model == "boltz2":
        return build_boltz_yaml(row, input_dir / f"{design_id}.yaml", input_field=input_field)
    if spec.parser == "summary" and model == "opendde":
        return build_opendde_fasta(row, input_dir / f"{design_id}.fasta", input_field=input_field)
    if spec.parser == "esmfold2":
        return build_esmfold2_fasta(row, input_dir / f"{design_id}.fasta", input_field=input_field)
    raise ValueError(f"no input builder for {model}")


def _parse_model(
    model_dir: Path,
    design_id: str,
    model: str,
    state: str,
    *,
    indexed_files: list[Path] | None = None,
) -> dict[str, object]:
    spec = get_predictor(model)
    if spec.parser == "esmfold2":
        return parse_esmfold(model_dir, design_id, prefix=_prefix("esmfold2", state), native=True)
    parsed = parse_summary_model(model_dir, design_id, model, indexed_files=indexed_files)
    prefix = f"{model}_"
    state_prefix = f"{_prefix(model, state)}_"
    return {
        (state_prefix + key[len(prefix):]) if key.startswith(prefix) else key: value
        for key, value in parsed.items()
    }


def _stage_empty(model: str, state: str, status: str, error: str) -> dict[str, object]:
    prefix = _prefix(model, state)
    return {
        f"{prefix}_status": status,
        f"{prefix}_error": error,
        f"{prefix}_log_path": "",
    }


def _run_one(
    row: dict[str, str],
    model: str,
    state: str,
    run_dir: Path,
    shard_dir: Path,
    config: dict[str, object],
    *,
    resume: bool,
) -> dict[str, object]:
    design_id = row["design_id"]
    model_root = shard_dir / "models" / model
    input_root = shard_dir / "inputs" / model
    model_dir = model_root / state / design_id if state != "single" else model_root / design_id
    input_dir = input_root / state if state != "single" else input_root
    model_dir.mkdir(parents=True, exist_ok=True)
    input_dir.mkdir(parents=True, exist_ok=True)
    try:
        input_path = _build_input(row, model, input_dir, state)
    except Exception as exc:
        return _stage_empty(model, state, "missing", f"input: {exc}")

    parsed = _parse_model(model_dir, design_id, model, state)
    prefix = _prefix(model, state)
    if resume and parsed.get(f"{prefix}_status") == "success":
        parsed[f"{prefix}_reused"] = "1"
        parsed[f"{prefix}_log_path"] = ""
        parsed[f"{prefix}_input_path"] = input_path
        return parsed

    log_name = f"level2_{model}_{state}_{design_id}.log" if state != "single" else f"level2_{model}_{design_id}.log"
    log_path = run_dir / "logs" / log_name
    context: dict[str, object] = {
        "run_dir": run_dir,
        "shard_dir": shard_dir,
        "design_id": design_id,
        "model": model,
        "state": state,
        "cid_state": state,
        "input": input_path,
        "outdir": model_dir,
        "model_input": input_path,
        "model_outdir": model_dir,
        "protenix_input": input_path if model == "protenix" else "",
        "protenix_input_dir": input_dir if model == "protenix" else "",
        "protenix_outdir": model_dir if model == "protenix" else "",
        "boltz_input": input_path if model == "boltz2" else "",
        "boltz_input_dir": input_dir if model == "boltz2" else "",
        "boltz_outdir": model_dir if model == "boltz2" else "",
        "opendde_input": input_path if model == "opendde" else "",
        "opendde_input_dir": input_dir if model == "opendde" else "",
        "opendde_outdir": model_dir if model == "opendde" else "",
        "esmfold2_input": input_path if model == "esmfold2" else "",
        "esmfold2_input_dir": input_dir if model == "esmfold2" else "",
        "esmfold2_outdir": model_dir if model == "esmfold2" else "",
    }
    command = _default_command(config, model, context)
    if command is None:
        return _stage_empty(model, state, "missing", f"{model} command is not configured")
    rc, tail = run_logged(command, log_path, cwd=shard_dir)
    parsed = _parse_model(model_dir, design_id, model, state)
    parsed[f"{prefix}_log_path"] = str(log_path)
    parsed[f"{prefix}_input_path"] = input_path
    if rc != 0:
        message = f"command returncode={rc}"
        if tail:
            message = f"{message}; {tail}"
        if parsed.get(f"{prefix}_status") == "success":
            parsed[f"{prefix}_status"] = "partial"
        elif parsed.get(f"{prefix}_status") == "missing":
            parsed[f"{prefix}_status"] = "failed"
        parsed[f"{prefix}_error"] = "; ".join(
            item for item in (str(parsed.get(f"{prefix}_error", "")), message) if item
        )
    return parsed


def _batch_context(
    run: Path,
    shard_dir: Path,
    model: str,
    state: str,
    input_dir: Path,
    output_dir: Path,
) -> dict[str, object]:
    """Build the command context for one model-wide batch invocation."""

    context: dict[str, object] = {
        "run_dir": run,
        "shard_dir": shard_dir,
        "model": model,
        "state": state,
        "cid_state": state,
        "input": input_dir,
        "outdir": output_dir,
        "model_input": input_dir,
        "model_outdir": output_dir,
        "protenix_input": input_dir if model == "protenix" else "",
        "protenix_input_dir": input_dir if model == "protenix" else "",
        "protenix_outdir": output_dir if model == "protenix" else "",
        "boltz_input": input_dir if model == "boltz2" else "",
        "boltz_input_dir": input_dir if model == "boltz2" else "",
        "boltz_outdir": output_dir if model == "boltz2" else "",
        "opendde_input": input_dir if model == "opendde" else "",
        "opendde_input_dir": input_dir if model == "opendde" else "",
        "opendde_outdir": output_dir if model == "opendde" else "",
        "esmfold2_input": input_dir if model == "esmfold2" else "",
        "esmfold2_input_dir": input_dir if model == "esmfold2" else "",
        "esmfold2_outdir": output_dir if model == "esmfold2" else "",
    }
    return context


def _run_model_batch(
    rows: list[dict[str, str]],
    model: str,
    state: str,
    run: Path,
    shard_dir: Path,
    config: dict[str, object],
    *,
    resume: bool,
) -> dict[str, dict[str, object]]:
    """Run one predictor once for every design in a shard.

    Predictors receive a directory, which prevents repeated model loading and
    keeps co-folding screening aligned with the batch input contract. Their
    output directory is parsed recursively per design after the single
    invocation completes.
    """

    model_root = shard_dir / "models" / model
    input_dir = shard_dir / "inputs" / model
    if state != "single":
        model_root = model_root / state
        input_dir = input_dir / state
    model_root.mkdir(parents=True, exist_ok=True)
    input_dir.mkdir(parents=True, exist_ok=True)
    # Index the batch output tree once.  The old per-design parser repeatedly
    # walked the complete nested predictor tree, which became quadratic for
    # large campaigns.
    indexed_files = [path for path in model_root.rglob("*") if path.is_file()]

    input_paths: dict[str, str] = {}
    result: dict[str, dict[str, object]] = {}
    for row in rows:
        design_id = row["design_id"]
        try:
            input_paths[design_id] = _build_input(row, model, input_dir, state)
        except Exception as exc:
            result[design_id] = _stage_empty(model, state, "missing", f"input: {exc}")

    pending = [row for row in rows if row["design_id"] not in result]
    prefix = _prefix(model, state)
    if not pending:
        return result

    # A complete batch can be reused without relaunching the model.  Partial
    # batches are relaunched so all designs share the same predictor settings.
    if resume:
        parsed_all = {
            row["design_id"]: _parse_model(
                model_root, row["design_id"], model, state, indexed_files=indexed_files
            )
            for row in pending
        }
        if all(parsed.get(f"{prefix}_status") == "success" for parsed in parsed_all.values()):
            for design_id, parsed in parsed_all.items():
                parsed[f"{prefix}_reused"] = "1"
                parsed[f"{prefix}_log_path"] = ""
                parsed[f"{prefix}_input_path"] = input_paths[design_id]
                result[design_id] = parsed
            return result

    log_name = (
        f"level2_{model}_{state}_batch.log"
        if state != "single"
        else f"level2_{model}_batch.log"
    )
    log_path = run / "logs" / log_name
    context = _batch_context(run, shard_dir, model, state, input_dir, model_root)
    command = _default_command(config, model, context)
    if command is None:
        for row in pending:
            result[row["design_id"]] = _stage_empty(
                model, state, "missing", f"{model} command is not configured"
            )
        return result
    rc, tail = run_logged(command, log_path, cwd=shard_dir)
    indexed_files = [path for path in model_root.rglob("*") if path.is_file()]
    for row in pending:
        design_id = row["design_id"]
        parsed = _parse_model(
            model_root, design_id, model, state, indexed_files=indexed_files
        )
        parsed[f"{prefix}_log_path"] = str(log_path)
        parsed[f"{prefix}_input_path"] = input_paths[design_id]
        if rc != 0:
            message = f"batch command returncode={rc}"
            if tail:
                message = f"{message}; {tail}"
            if parsed.get(f"{prefix}_status") == "success":
                parsed[f"{prefix}_status"] = "partial"
            elif parsed.get(f"{prefix}_status") == "missing":
                parsed[f"{prefix}_status"] = "failed"
            parsed[f"{prefix}_error"] = "; ".join(
                item
                for item in (str(parsed.get(f"{prefix}_error", "")), message)
                if item
            )
        result[design_id] = parsed
    return result


def _cid_enabled(config: dict[str, object], explicit: bool | None) -> bool:
    if explicit is not None:
        return explicit
    mode = str(config.get("level2_mode", config.get("cid_mode", "single"))).strip().lower()
    return mode in {"cid", "apo_holo", "holo_apo"}


def run_shard(
    run_dir: str | Path,
    shard_index: int,
    *,
    resume: bool = True,
    cid: bool | None = None,
) -> list[dict[str, object]]:
    run = Path(run_dir).resolve()
    config = load_config(run)
    rows = load_shard(run, shard_index)
    shard_dir = run / "artifacts" / f"shard_{shard_index:05d}"
    shard_dir.mkdir(parents=True, exist_ok=True)
    states = CID_STATES if _cid_enabled(config, cid) else ("single",)
    models = configured_models(config)
    output: list[dict[str, object]] = []
    batch_results: dict[tuple[str, str], dict[str, dict[str, object]]] = {}
    timings: list[dict[str, object]] = []
    for model in models:
        for state in states:
            started = time.perf_counter()
            batch_results[(model, state)] = _run_model_batch(
                rows, model, state, run, shard_dir, config, resume=resume
            )
            finished = time.perf_counter()
            values = batch_results[(model, state)]
            statuses = [str(item.get(f"{_prefix(model, state)}_status", "missing")) for item in values.values()]
            timings.append({
                "shard_index": shard_index,
                "model": model,
                "state": state,
                "design_count": len(rows),
                "result_count": len(values),
                "success_count": sum(status == "success" for status in statuses),
                "elapsed_seconds": round(finished - started, 3),
                "resumed": bool(resume and values and all(
                    item.get(f"{_prefix(model, state)}_reused") == "1" for item in values.values()
                )),
            })
    for row in rows:
        result: dict[str, object] = {"design_id": row["design_id"]}
        for model in models:
            for state in states:
                result.update(batch_results[(model, state)].get(row["design_id"], {}))
        statuses = {
            _prefix(model, state): str(result.get(f"{_prefix(model, state)}_status", "missing"))
            for model in models for state in states
        }
        result["status"] = status_from_models(statuses)
        result["error"] = "; ".join(
            str(result.get(f"{_prefix(model, state)}_error", ""))
            for model in models for state in states
            if result.get(f"{_prefix(model, state)}_error", "")
        )
        result["task_index"] = row.get("task_index", "")
        result["shard_index"] = str(shard_index)
        output.append(result)
    # Generate a shard-local monitor report after all model outputs and parser
    # states exist.  Keeping it shard-local avoids races between Slurm array
    # workers; aggregate_level2.py later writes the run-level report.
    for state in states:
        monitor_root = shard_dir / "models" / "esmfold2"
        if state != "single":
            monitor_root = monitor_root / state
        esmfold2_rows = scan_esmfold2(monitor_root)
        monitor_name = (
            f"esmfold2_monitor_{state}_shard_{shard_index:05d}.csv"
            if state != "single"
            else f"esmfold2_monitor_shard_{shard_index:05d}.csv"
        )
        write_esmfold2_csv(
            esmfold2_rows,
            run / "metrics" / monitor_name,
        )
    write_part(run, shard_index, output)
    # Keep timing data shard-local so concurrent Slurm workers never contend
    # over one file.  The aggregate stage creates the run-level table.
    atomic_json(shard_dir / "stage_timings.json", timings)
    return output


def main() -> int:
    parser = argparse.ArgumentParser(description="Run one ProFold Level 2 shard")
    parser.add_argument("--run-dir", required=True)
    parser.add_argument("--shard-index", type=int, required=True)
    parser.add_argument("--no-resume", action="store_true")
    parser.add_argument("--cid", action="store_true", help="run apo and holo states")
    parser.add_argument("--no-cid", action="store_true", help="force single-state mode")
    args = parser.parse_args()
    if args.cid and args.no_cid:
        parser.error("--cid and --no-cid are mutually exclusive")
    explicit_cid = True if args.cid else False if args.no_cid else None
    run_shard(args.run_dir, args.shard_index, resume=not args.no_resume, cid=explicit_cid)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
