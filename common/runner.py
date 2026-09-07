#!/usr/bin/env python3
"""Command, logging and run-directory primitives used by all levels."""

from __future__ import annotations

import json
import shlex
import subprocess
from pathlib import Path
from string import Formatter
from typing import Mapping, Sequence


def atomic_json(path: str | Path, payload: object) -> None:
    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    tmp = output.with_name(f".{output.name}.tmp")
    tmp.write_text(json.dumps(payload, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    tmp.replace(output)


def read_json(path: str | Path, default: object | None = None) -> object:
    input_path = Path(path)
    if not input_path.is_file():
        if default is not None:
            return default
        raise FileNotFoundError(input_path)
    return json.loads(input_path.read_text(encoding="utf-8"))


def tail_text(path: str | Path, lines: int = 20) -> str:
    input_path = Path(path)
    if not input_path.is_file():
        return ""
    content = input_path.read_text(encoding="utf-8", errors="replace").splitlines()
    return "\n".join(content[-lines:])


def command_template(value: object, context: Mapping[str, object]) -> list[str] | None:
    if value in (None, ""):
        return None
    tokens = shlex.split(value) if isinstance(value, str) else [str(item) for item in value]  # type: ignore[arg-type]
    formatter = Formatter()
    rendered: list[str] = []
    for token in tokens:
        fields = [name for _, name, _, _ in formatter.parse(token) if name]
        missing = [field for field in fields if field not in context]
        if missing:
            raise KeyError(f"command placeholder(s) missing: {', '.join(missing)}")
        rendered.append(token.format_map({key: str(val) for key, val in context.items()}))
    return rendered


def configured_command(config: Mapping[str, object], name: str, context: Mapping[str, object]) -> list[str] | None:
    commands = config.get("commands", {})
    if isinstance(commands, Mapping) and name in commands:
        return command_template(commands[name], context)
    return None


def run_logged(
    command: Sequence[str],
    log_path: str | Path,
    *,
    cwd: str | Path | None = None,
    env: Mapping[str, str] | None = None,
) -> tuple[int, str]:
    log = Path(log_path)
    log.parent.mkdir(parents=True, exist_ok=True)
    with log.open("w", encoding="utf-8") as handle:
        handle.write("command=" + shlex.join([str(item) for item in command]) + "\n")
        handle.flush()
        try:
            completed = subprocess.run(
                [str(item) for item in command],
                cwd=str(cwd) if cwd else None,
                env=dict(env) if env else None,
                stdout=handle,
                stderr=subprocess.STDOUT,
                check=False,
            )
            handle.write(f"\nreturncode={completed.returncode}\n")
            return completed.returncode, tail_text(log)
        except OSError as exc:
            handle.write(f"\nlauncher_error={exc!r}\n")
            return 127, tail_text(log)


def stage_result(model: str, status: str, *, error: str = "", log_path: str = "", **extra: object) -> dict[str, object]:
    result: dict[str, object] = {
        f"{model}_status": status,
        f"{model}_error": error,
        f"{model}_log_path": log_path,
    }
    result.update(extra)
    return result

