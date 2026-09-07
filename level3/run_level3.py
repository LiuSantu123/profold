#!/usr/bin/env python3
"""Run one Level 3 shard with the historical NetSolP and TemBERTure scripts."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

V37 = Path(__file__).resolve().parents[1]
if str(V37) not in sys.path:
    sys.path.insert(0, str(V37))

from common.input_builder import protein_records, status_from_models  # noqa: E402
from common.manifest import metadata  # noqa: E402
from common.model_parsers import parse_netsolp, parse_temberture  # noqa: E402
from common.pipeline import load_config, load_shard, write_part  # noqa: E402
from common.runner import command_template, configured_command, run_logged  # noqa: E402


NETSOLP_WRAPPER = V37 / "wrappers" / "run_netsolp.sh"
TEMBERTURE_SCRIPT = Path(
    "/xcfhome/yhliu/14_magpcr/06_zqq_5cpm/7_tm_pre/1_temberture/run_temberture_v7.py"
)


def _analysis_sequence(row: dict[str, str], config: dict[str, object]) -> str:
    """Resolve the sequence scored by sequence-only models.

    Single-chain rows are used directly.  For a multichain row, an explicit
    ``level3_sequence``/metadata value is preferred; otherwise the config must
    opt into ``first`` or ``concat`` to avoid silently scoring the wrong chain.
    """

    meta = metadata(row)
    explicit = row.get("level3_sequence", "").strip() or str(
        meta.get("level3_sequence", meta.get("sequence_for_level3", ""))
    ).strip()
    if explicit:
        return explicit.replace(" ", "").replace("\n", "").upper()
    records = protein_records(row)
    if len(records) == 1:
        return records[0][1]
    mode = str(config.get("level3_multichain_mode", "error")).lower()
    if mode == "first":
        return records[0][1]
    if mode == "concat":
        return "".join(sequence for _chain, sequence in records)
    raise ValueError(
        "multiple protein chains; provide level3_sequence or set "
        "level3_multichain_mode to first/concat"
    )


def _write_fasta(rows: list[dict[str, str]], sequences: dict[str, str], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(f".{path.name}.tmp")
    with tmp.open("w", encoding="utf-8") as handle:
        for row in rows:
            design_id = row["design_id"]
            sequence = sequences.get(design_id, "")
            if sequence:
                handle.write(f">{design_id}\n{sequence}\n")
    tmp.replace(path)


def _stage_command(
    config: dict[str, object], name: str, context: dict[str, object]
) -> list[str] | None:
    configured = configured_command(config, name, context)
    if configured is not None:
        return configured
    if name == "netsolp":
        wrapper = str(config.get("netsolp_wrapper", NETSOLP_WRAPPER))
        return command_template(
            ["bash", wrapper, "{level3_fasta}", "{netsolp_raw}", "{netsolp_csv}"],
            context,
        )
    if name == "temberture":
        python_bin = str(config.get("temberture_python", sys.executable))
        script = str(config.get("temberture_script", TEMBERTURE_SCRIPT))
        return command_template(
            [python_bin, script, "{level3_fasta}", "--output", "{temberture_csv}"],
            context,
        )
    raise ValueError(name)


def _append_command_error(
    parsed: dict[str, object], model: str, rc: int, tail: str
) -> dict[str, object]:
    if rc == 0:
        return parsed
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


def _run_stage(
    model: str,
    rows: list[dict[str, str]],
    sequences: dict[str, str],
    run: Path,
    shard_dir: Path,
    config: dict[str, object],
    *,
    resume: bool,
) -> dict[str, dict[str, object]]:
    output_dir = shard_dir / "level3"
    output_dir.mkdir(parents=True, exist_ok=True)
    fasta_path = output_dir / "inputs.fasta"
    _write_fasta(rows, sequences, fasta_path)
    if model == "netsolp":
        csv_path = output_dir / "netsolp_extracted.csv"
        raw_path = output_dir / "netsolp_raw.csv"
        parser = parse_netsolp
    else:
        csv_path = output_dir / "temberture.csv"
        raw_path = output_dir / "temberture.log"
        parser = parse_temberture

    valid_rows = [row for row in rows if row["design_id"] in sequences]
    parsed = {row["design_id"]: parser(csv_path, row["design_id"]) for row in valid_rows}
    need_run = not resume or any(
        item.get(f"{model}_status") != "success" for item in parsed.values()
    )
    if need_run and valid_rows:
        log_path = run / "logs" / f"level3_{model}_shard_{shard_dir.name.split('_')[-1]}.log"
        context: dict[str, object] = {
            "run_dir": run,
            "shard_dir": shard_dir,
            "level3_fasta": fasta_path,
            "fasta": fasta_path,
            "input": fasta_path,
            "netsolp_raw": raw_path,
            "netsolp_csv": csv_path,
            "temberture_csv": csv_path,
            "output": csv_path,
            "outdir": output_dir,
        }
        command = _stage_command(config, model, context)
        if command is not None:
            rc, tail = run_logged(command, log_path, cwd=shard_dir)
            parsed = {row["design_id"]: parser(csv_path, row["design_id"]) for row in valid_rows}
            for item in parsed.values():
                item[f"{model}_log_path"] = str(log_path)
                item = _append_command_error(item, model, rc, tail)
        else:
            parsed = {
                row["design_id"]: {
                    f"{model}_status": "missing",
                    f"{model}_error": f"{model} command is not configured",
                    f"{model}_log_path": "",
                }
                for row in valid_rows
            }
    for item in parsed.values():
        item[f"{model}_input_fasta"] = str(fasta_path)
    return parsed


def run_shard(run_dir: str | Path, shard_index: int, *, resume: bool = True) -> list[dict[str, object]]:
    run = Path(run_dir).resolve()
    config = load_config(run)
    rows = load_shard(run, shard_index)
    shard_dir = run / "artifacts" / f"shard_{shard_index:05d}"
    shard_dir.mkdir(parents=True, exist_ok=True)
    sequences: dict[str, str] = {}
    sequence_errors: dict[str, str] = {}
    for row in rows:
        try:
            sequences[row["design_id"]] = _analysis_sequence(row, config)
        except Exception as exc:
            sequence_errors[row["design_id"]] = str(exc)

    models = ("netsolp", "temberture")
    stage_results = {
        model: _run_stage(model, rows, sequences, run, shard_dir, config, resume=resume)
        for model in models
    }
    output: list[dict[str, object]] = []
    for row in rows:
        design_id = row["design_id"]
        result: dict[str, object] = {"design_id": design_id}
        if design_id in sequence_errors:
            for model in models:
                result.update(
                    {
                        f"{model}_status": "missing",
                        f"{model}_error": f"level3 sequence: {sequence_errors[design_id]}",
                        f"{model}_log_path": "",
                    }
                )
        else:
            for model in models:
                result.update(stage_results[model].get(design_id, {}))
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
    parser = argparse.ArgumentParser(description="Run one v37 Level 3 shard")
    parser.add_argument("--run-dir", required=True)
    parser.add_argument("--shard-index", type=int, required=True)
    parser.add_argument("--no-resume", action="store_true")
    args = parser.parse_args()
    run_shard(args.run_dir, args.shard_index, resume=not args.no_resume)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
