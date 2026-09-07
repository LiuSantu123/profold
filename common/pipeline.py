#!/usr/bin/env python3
"""Run preparation, shard handling and aggregation shared by the three levels."""

from __future__ import annotations

import math
from pathlib import Path
from typing import Iterable, Mapping

from .manifest import read_manifest, row_fingerprint, write_tsv
from .runner import atomic_json, read_json


RUN_DIRS = ("input", "tasks", "results/parts", "metrics", "status/parts", "artifacts", "logs")


def prepare_run(
    input_manifest: str | Path,
    outdir: str | Path,
    *,
    level: str,
    shard_size: int,
    config: Mapping[str, object] | None = None,
) -> tuple[Path, list[dict[str, str]]]:
    if shard_size < 1:
        raise ValueError("shard_size must be >= 1")
    input_path = Path(input_manifest).expanduser().resolve()
    output = Path(outdir).expanduser().resolve()
    rows = read_manifest(input_path, require_sequence=True)
    if level == 'level1' and str((config or {}).get('af3_mode', '')).strip().lower() == 'cid':
        from level1.cid import validate_cid_rows
        validate_cid_rows(rows)
    output.mkdir(parents=True, exist_ok=True)
    for relative in RUN_DIRS:
        (output / relative).mkdir(parents=True, exist_ok=True)

    normalized = []
    for index, row in enumerate(rows):
        item = dict(row)
        item["task_index"] = str(index)
        item["shard_index"] = str(index // shard_size)
        item["input_fingerprint"] = row_fingerprint(row)
        normalized.append(item)

    write_tsv(output / "input/input_manifest.tsv", rows)
    write_tsv(output / "tasks/task_manifest.tsv", normalized, preferred=("task_index", "shard_index"))
    n_shards = math.ceil(len(normalized) / shard_size)
    for shard_index in range(n_shards):
        shard_rows = [row for row in normalized if int(row["shard_index"]) == shard_index]
        write_tsv(output / f"tasks/shard_{shard_index:05d}.tsv", shard_rows, preferred=("task_index", "shard_index"))
    (output / "tasks/n_shards.txt").write_text(f"{n_shards}\n", encoding="utf-8")
    payload = dict(config or {})
    payload.update(
        {
            "level": level,
            "input_manifest": str(input_path),
            "outdir": str(output),
            "shard_size": shard_size,
            "n_tasks": len(rows),
            "n_shards": n_shards,
        }
    )
    atomic_json(output / "config.json", payload)
    return output, normalized


def load_config(run_dir: str | Path) -> dict[str, object]:
    value = read_json(Path(run_dir) / "config.json", default={})
    if not isinstance(value, dict):
        raise ValueError(f"run config must be an object: {Path(run_dir) / 'config.json'}")
    return dict(value)


def load_shard(run_dir: str | Path, shard_index: int) -> list[dict[str, str]]:
    from .manifest import read_manifest

    run = Path(run_dir).resolve()
    return read_manifest(run / f"tasks/shard_{shard_index:05d}.tsv", require_sequence=True)


def part_path(run_dir: str | Path, shard_index: int) -> Path:
    return Path(run_dir).resolve() / f"results/parts/part_{shard_index:05d}.tsv"


def write_part(run_dir: str | Path, shard_index: int, rows: Iterable[Mapping[str, object]]) -> Path:
    output = part_path(run_dir, shard_index)
    write_tsv(output, rows, preferred=("design_id", "status", "error"))
    return output


def read_part_files(run_dir: str | Path) -> list[Path]:
    return sorted((Path(run_dir).resolve() / "results/parts").glob("part_*.tsv"))


def join_results(
    run_dir: str | Path,
    *,
    model_names: Iterable[str],
    missing_error: str = "worker part missing",
) -> list[dict[str, str]]:
    from .manifest import read_manifest

    run = Path(run_dir).resolve()
    base_rows = read_manifest(run / "tasks/task_manifest.tsv", require_sequence=True)
    by_id: dict[str, dict[str, str]] = {}
    for part in read_part_files(run):
        for row in read_manifest(part):
            design_id = row["design_id"]
            if design_id in by_id:
                raise ValueError(f"duplicate result for design_id={design_id}")
            by_id[design_id] = row
    output: list[dict[str, str]] = []
    for base in base_rows:
        row = dict(base)
        result = by_id.get(base["design_id"])
        if result is None:
            row.update({"status": "missing", "error": missing_error})
            for model in model_names:
                row[f"{model}_status"] = "missing"
                row[f"{model}_error"] = missing_error
        else:
            # Worker parts intentionally contain only result fields.  They are
            # read through the normal manifest loader, which adds blank core
            # fields (sequence, input paths, ...).  Do not let those defaults
            # erase the original input row when building a downstream manifest.
            row.update(
                {
                    key: value
                    for key, value in result.items()
                    if key not in base or value != ""
                }
            )
        output.append(row)
    return output


def status_rows(rows: Iterable[Mapping[str, object]], level: str, model_names: Iterable[str]) -> list[dict[str, object]]:
    models = list(model_names)
    output: list[dict[str, object]] = []
    for row in rows:
        base = {"design_id": row.get("design_id", ""), "level": level, "status": row.get("status", ""), "error": row.get("error", "")}
        output.append(base)
        for model in models:
            output.append(
                {
                    "design_id": row.get("design_id", ""),
                    "level": level,
                    "model": model,
                    "status": row.get(f"{model}_status", "missing"),
                    "error": row.get(f"{model}_error", ""),
                    "log_path": row.get(f"{model}_log_path", ""),
                }
            )
    return output


def metric_rows(rows: Iterable[Mapping[str, object]], model_names: Iterable[str]) -> list[dict[str, object]]:
    """Emit one auditable row per design/model without filling missing scores."""

    models = list(model_names)
    result: list[dict[str, object]] = []
    for row in rows:
        for model in models:
            prefix = f"{model}_"
            metrics = {
                key[len(prefix) :]: value
                for key, value in row.items()
                if key.startswith(prefix)
                and not key.endswith(("_status", "_error", "_log_path"))
                and key not in {f"{model}_output_dir"}
            }
            item: dict[str, object] = {
                "design_id": row.get("design_id", ""),
                "model": model,
                "status": row.get(f"{model}_status", "missing"),
            }
            item.update(metrics)
            result.append(item)
    return result
