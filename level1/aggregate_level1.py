#!/usr/bin/env python3
"""Merge Level 1 shard results into a complete result/status manifest."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

V37 = Path(__file__).resolve().parents[1]
if str(V37) not in sys.path:
    sys.path.insert(0, str(V37))

from common.manifest import write_csv, write_tsv  # noqa: E402
from common.pipeline import join_results, metric_rows, status_rows, load_config  # noqa: E402


def main() -> int:
    parser = argparse.ArgumentParser(description="Aggregate v37 Level 1")
    parser.add_argument("--run-dir", required=True)
    args = parser.parse_args()
    run = Path(args.run_dir).resolve()
    rows = join_results(run, model_names=("esmfold", "af3"))
    write_tsv(run / "results/result_manifest.tsv", rows)
    write_csv(run / "metrics/level1.csv", metric_rows(rows, ("esmfold", "af3")))
    models = ["esmfold", "af3"]
    if load_config(run).get('af3_mode') in {'cid', 'for_wj_cid'}:
        for state in ('apo', 'holo'):
            model = f'af3_{state}'
            models.append(model)
            write_csv(run / f'metrics/cid_{state}.csv', metric_rows(rows, (model,)))
    write_tsv(run / "status/status.tsv", status_rows(rows, "level1", models))
    counts: dict[str, int] = {}
    for row in rows:
        counts[str(row.get("status", "missing"))] = counts.get(str(row.get("status", "missing")), 0) + 1
    print(f"aggregated level1: {len(rows)} rows -> {run / 'results/result_manifest.tsv'}; {counts}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
