#!/usr/bin/env python3
"""Run one Level 1 shard: batch ESMFold, then batch AF3."""

from __future__ import annotations

import argparse
import hashlib
import sys
from pathlib import Path

V37 = Path(__file__).resolve().parents[1]
if str(V37) not in sys.path:
    sys.path.insert(0, str(V37))

from common.input_builder import build_af3_json, protein_records, status_from_models  # noqa: E402
from common.manifest import write_fasta  # noqa: E402
from common.model_parsers import parse_af3, parse_esmfold  # noqa: E402
from common.pipeline import load_config, load_shard, write_part  # noqa: E402
from common.runner import configured_command, run_logged  # noqa: E402
from level1.cid import cid_command, validate_cid_rows  # noqa: E402


ESMFOLD_BIN = Path("/xcfhome/yzmeng/miniconda3/envs/zb/bin/esm-fold")
AF3_MODES = {"direct", "for_wj_cid", "cid"}


def _seed(design_id: str) -> int:
    return int(hashlib.sha256(design_id.encode("utf-8")).hexdigest()[:8], 16) % 100000 + 1


def _af3_mode(config: dict[str, object]) -> str:
    mode = str(config.get("af3_mode", "direct")).strip().lower()
    if mode not in AF3_MODES:
        choices = ", ".join(sorted(AF3_MODES))
        raise ValueError(f"unsupported af3_mode={mode!r}; choose one of: {choices}")
    return mode


def _validate_for_wj_cid_row(row: dict[str, str]) -> None:
    records = protein_records(row)
    if not records:
        raise ValueError(
            "for_wj_cid requires at least one protein chain; "
            f"design {row['design_id']} has {len(records)}"
        )


def _write_af3_fasta(rows: list[dict[str, str]], path: Path, mode: str) -> int:
    if mode not in {"for_wj_cid", "cid"}:
        return write_fasta(rows, path)
    path.parent.mkdir(parents=True, exist_ok=True)
    count = 0
    chain_order = None
    tmp = path.with_name(f".{path.name}.tmp")
    with tmp.open("w", encoding="utf-8") as handle:
        for row in rows:
            records = protein_records(row)
            ids = [cid for cid, _ in records]
            if chain_order is not None and ids != chain_order:
                raise ValueError("for_wj_cid shard requires one shared protein chain order/template")
            chain_order = ids
            if not records:
                raise ValueError(
                    "for_wj_cid requires at least one protein chain; "
                    f"design {row['design_id']} has {len(records)}"
                )
            handle.write(f">{row['design_id']}\n{':'.join(seq for _, seq in records)}\n")
            count += 1
    tmp.replace(path)
    return count


def _default_command(config: dict[str, object], name: str, context: dict[str, object]) -> list[str] | None:
    configured = configured_command(config, name, context)
    if configured is not None:
        return configured
    if name == "af3" and context.get("af3_mode") == "cid":
        return cid_command(config)
    if name == "esmfold":
        command = [
            str(config.get("esmfold_bin", ESMFOLD_BIN)),
            "-i", "{esm_input}", "-o", "{esm_outdir}",
        ]
        for key, flag in (
            ("esmfold_model_dir", "-m"),
            ("esmfold_num_recycles", "--num-recycles"),
            ("esmfold_max_tokens_per_batch", "--max-tokens-per-batch"),
            ("esmfold_chunk_size", "--chunk-size"),
        ):
            if config.get(key) is not None:
                command.extend([flag, str(config[key])])
        for key, flag in (("esmfold_cpu_only", "--cpu-only"), ("esmfold_cpu_offload", "--cpu-offload")):
            if config.get(key, False):
                command.append(flag)
        return command
    # The for_wj wrapper owns AF3 JSON generation and must be explicitly wired
    # through commands.af3; its CLI is not interchangeable with run_alphafold.
    if name == "af3" and context.get("af3_mode") == "for_wj_cid":
        return None

    # AF3 model/database locations differ between installations.  Require
    # them in config rather than launching a command with an invalid path.
    script = config.get("af3_script")
    model_dir = config.get("af3_model_dir")
    db_dir = config.get("af3_db_dir")
    if script and model_dir and db_dir:
        python_bin = str(config.get("af3_python", sys.executable))
        return [
            python_bin,
            str(script),
            "--input_dir",
            "{af3_input_dir}",
            "--output_dir",
            "{af3_outdir}",
            "--model_dir",
            "{af3_model_dir}",
            "--db_dir",
            "{af3_db_dir}",
        ]
    return None


def _context(
    run_dir: Path,
    shard_dir: Path,
    rows: list[dict[str, str]],
    config: dict[str, object],
) -> dict[str, object]:
    af3_dir = shard_dir / "af3"
    return {
        "run_dir": run_dir,
        "shard_dir": shard_dir,
        "esm_input": shard_dir / "inputs.fasta",
        "esm_outdir": shard_dir / "esmfold",
        "af3_input_dir": shard_dir / "af3_inputs",
        "af3_outdir": af3_dir,
        "af3_fasta": shard_dir / "af3_inputs.fasta",
        "af3_csv": af3_dir / "af3_summary.csv",
        "af3_archive_dir": af3_dir / "archive",
        "af3_json_archive_dir": af3_dir / "json_archive",
        "af3_esm_outdir": shard_dir / "esmfold",
        "af3_sequence_chain_ids": [cid for cid, _ in protein_records(rows[0])] if rows else [],
        "cid_apo_json_path": rows[0].get("cid_apo_json_path", "") if rows else "",
        "cid_holo_json_path": rows[0].get("cid_holo_json_path", "") if rows else "",
        "af3_temp_dir": shard_dir / ".af3_temp",
        "af3_model_dir": "",
        "af3_db_dir": "",
        "esmfold_bin": config.get("esmfold_bin", ESMFOLD_BIN),
        "n_tasks": len(rows),
    }


def _run_stage(
    config: dict[str, object],
    name: str,
    context: dict[str, object],
    log_path: Path,
) -> tuple[int | None, str]:
    command = _default_command(config, name, context)
    if command is None:
        return None, f"{name} command is not configured"
    rendered = configured_command(config, name, context)
    if rendered is None:
        rendered = [str(item).format_map({key: str(value) for key, value in context.items()}) for item in command]
    if name == "af3" and context.get("af3_mode") in {"for_wj_cid", "cid"} and "--sequence-chain-ids" not in rendered:
        rendered.extend(["--sequence-chain-ids", *context["af3_sequence_chain_ids"]])
    rc, tail = run_logged(rendered, log_path, cwd=context["shard_dir"])
    return rc, tail


def _empty_stage(model: str, status: str, error: str) -> dict[str, object]:
    return {
        f"{model}_status": status,
        f"{model}_error": error,
        f"{model}_log_path": "",
    }


def run_shard(run_dir: str | Path, shard_index: int, *, resume: bool = True) -> list[dict[str, object]]:
    run = Path(run_dir).resolve()
    config = load_config(run)
    af3_mode = _af3_mode(config)
    rows = load_shard(run, shard_index)
    if af3_mode == "cid":
        validate_cid_rows(rows)
    shard_dir = run / "artifacts" / f"shard_{shard_index:05d}"
    input_dir = shard_dir / "af3_inputs"
    esm_dir = shard_dir / "esmfold"
    af3_dir = shard_dir / "af3"
    input_dir.mkdir(parents=True, exist_ok=True)
    shard_dir.mkdir(parents=True, exist_ok=True)

    results: dict[str, dict[str, object]] = {row["design_id"]: {"design_id": row["design_id"]} for row in rows}
    valid_af3: list[dict[str, str]] = []
    for row in rows:
        design_id = row["design_id"]
        try:
            if af3_mode == "direct":
                build_af3_json(row, input_dir / f"{design_id}.json", seed=_seed(design_id))
            else:
                _validate_for_wj_cid_row(row)
            valid_af3.append(row)
        except Exception as exc:
            results[design_id].update(_empty_stage("af3", "missing", f"AF3 input: {exc}"))
    esm_rows = [dict(row, sequence=":".join(seq for _, seq in protein_records(row))) for row in rows]
    write_fasta(esm_rows, shard_dir / "inputs.fasta")
    _write_af3_fasta(valid_af3, shard_dir / "af3_inputs.fasta", af3_mode)

    context = _context(run, shard_dir, rows, config)
    context["af3_mode"] = af3_mode
    context["af3_model_dir"] = config.get("af3_model_dir", "")
    context["af3_db_dir"] = config.get("af3_db_dir", "")
    esm_log = run / "logs" / f"level1_esmfold_shard_{shard_index:05d}.log"
    af3_log = run / "logs" / f"level1_af3_shard_{shard_index:05d}.log"

    native_esm = configured_command(config, "esmfold", context) is None
    esm_existing = {row["design_id"]: parse_esmfold(esm_dir, row["design_id"], native=native_esm) for row in rows}
    need_esm = not resume or any(value["esmfold_status"] != "success" for value in esm_existing.values())
    esm_rc: int | None = 0
    esm_launch_error = ""
    if need_esm:
        try:
            esm_rc, esm_launch_error = _run_stage(config, "esmfold", context, esm_log)
        except Exception as exc:
            esm_rc, esm_launch_error = 127, repr(exc)
    for row in rows:
        design_id = row["design_id"]
        parsed = parse_esmfold(esm_dir, design_id, native=native_esm)
        if not need_esm and parsed["esmfold_status"] == "success":
            parsed["esmfold_reused"] = "1"
        if need_esm and esm_rc not in (None, 0) and parsed["esmfold_status"] == "success":
            parsed["esmfold_status"] = "partial"
            parsed["esmfold_error"] = f"command returncode={esm_rc}; {esm_launch_error}".strip()
        elif need_esm and esm_rc is None:
            parsed["esmfold_status"] = "missing"
            parsed["esmfold_error"] = esm_launch_error
        parsed["esmfold_log_path"] = str(esm_log) if need_esm else ""
        results[design_id].update(parsed)

    af3_existing = {
        row["design_id"]: parse_af3(af3_dir, row["design_id"], mode=af3_mode)
        for row in valid_af3
    }
    need_af3 = bool(valid_af3) and (not resume or any(value["af3_status"] != "success" for value in af3_existing.values()))
    af3_rc: int | None = 0
    af3_launch_error = ""
    if need_af3:
        try:
            af3_rc, af3_launch_error = _run_stage(config, "af3", context, af3_log)
        except Exception as exc:
            af3_rc, af3_launch_error = 127, repr(exc)
    for row in valid_af3:
        design_id = row["design_id"]
        parsed = parse_af3(af3_dir, design_id, mode=af3_mode)
        if not need_af3 and parsed["af3_status"] == "success":
            parsed["af3_reused"] = "1"
        if need_af3 and af3_rc not in (None, 0) and parsed["af3_status"] == "success":
            parsed["af3_status"] = "partial"
            parsed["af3_error"] = f"command returncode={af3_rc}; {af3_launch_error}".strip()
        elif need_af3 and af3_rc is None:
            parsed["af3_status"] = "missing"
            parsed["af3_error"] = af3_launch_error
        parsed["af3_log_path"] = str(af3_log) if need_af3 else ""
        results[design_id].update(parsed)

    output: list[dict[str, object]] = []
    for row in rows:
        result = results[row["design_id"]]
        statuses = {model: str(result.get(f"{model}_status", "missing")) for model in ("esmfold", "af3")}
        result["status"] = status_from_models(statuses)
        errors = [str(result.get(f"{model}_error", "")) for model in statuses if result.get(f"{model}_error", "")]
        result["error"] = "; ".join(errors)
        result["task_index"] = row.get("task_index", "")
        result["shard_index"] = str(shard_index)
        output.append(result)
    write_part(run, shard_index, output)
    return output


def main() -> int:
    parser = argparse.ArgumentParser(description="Run one v37 Level 1 shard")
    parser.add_argument("--run-dir", required=True)
    parser.add_argument("--shard-index", type=int, required=True)
    parser.add_argument("--no-resume", action="store_true")
    args = parser.parse_args()
    run_shard(args.run_dir, args.shard_index, resume=not args.no_resume)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
