#!/usr/bin/env python3
"""Prepare a Level 1 run directory and shard manifest."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

V37 = Path(__file__).resolve().parents[1]
if str(V37) not in sys.path:
    sys.path.insert(0, str(V37))

from common.pipeline import prepare_run  # noqa: E402


def main() -> int:
    parser = argparse.ArgumentParser(description="Prepare v37 Level 1 ESMFold + AF3 run")
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--shard-size", type=int, default=32)
    parser.add_argument("--config", help="JSON with commands and model paths")
    args = parser.parse_args()
    config = {}
    if args.config:
        config = json.loads(Path(args.config).read_text(encoding="utf-8"))
        if not isinstance(config, dict):
            raise SystemExit("--config must contain a JSON object")
    output, rows = prepare_run(
        args.manifest,
        args.outdir,
        level="level1",
        shard_size=args.shard_size,
        config=config,
    )
    print(f"prepared level1: {len(rows)} tasks, {output / 'tasks/n_shards.txt'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
