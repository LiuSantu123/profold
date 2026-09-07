#!/usr/bin/env python3
"""Aggregate Level 3 shard outputs."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

V37 = Path(__file__).resolve().parents[1]
if str(V37) not in sys.path:
    sys.path.insert(0, str(V37))

from common.manifest import write_csv, write_tsv  # noqa: E402
from common.pipeline import join_results, metric_rows, status_rows  # noqa: E402


def main() -> int:
    parser = argparse.ArgumentParser(description="Aggregate v37 Level 3")
    parser.add_argument("--run-dir", required=True)
    args = parser.parse_args()
    run = Path(args.run_dir).resolve()
    models = ("netsolp", "temberture")
    rows = join_results(run, model_names=models)
    write_tsv(run / "results/result_manifest.tsv", rows)
    write_csv(run / "metrics/level3.csv", metric_rows(rows, models))
    write_tsv(run / "status/status.tsv", status_rows(rows, "level3", models))
    counts: dict[str, int] = {}
    for row in rows:
        counts[str(row.get("status", "missing"))] = counts.get(str(row.get("status", "missing")), 0) + 1
    print(f"aggregated level3: {len(rows)} rows -> {run / 'results/result_manifest.tsv'}; {counts}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
