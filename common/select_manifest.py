#!/usr/bin/env python3
"""Explicitly select rows for the next level; never silently filter results."""

from __future__ import annotations

import argparse
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

V37 = Path(__file__).resolve().parents[1]
if str(V37) not in sys.path:
    sys.path.insert(0, str(V37))

from common.manifest import read_manifest, write_tsv  # noqa: E402


def _id_set(path: str | None, inline: str | None) -> list[str] | None:
    values: list[str] = []
    if path:
        values.extend(line.strip() for line in Path(path).read_text(encoding="utf-8").splitlines())
    if inline:
        values.extend(item.strip() for item in inline.split(","))
    values = [item for item in values if item and not item.startswith("#")]
    return values or None


def main() -> int:
    parser = argparse.ArgumentParser(description="Create a user-selected downstream manifest")
    parser.add_argument("--result-manifest", required=True)
    parser.add_argument("--ids-file")
    parser.add_argument("--ids", help="comma-separated design_id values")
    parser.add_argument("--top-n", type=int)
    parser.add_argument("--sort-field")
    parser.add_argument("--ascending", action="store_true")
    parser.add_argument("--require-success", action="store_true")
    parser.add_argument("-o", "--output", required=True)
    args = parser.parse_args()

    if bool(args.ids_file or args.ids) == bool(args.top_n is not None):
        parser.error("choose exactly one of --ids/--ids-file or --top-n")
    if args.top_n is not None and args.top_n < 1:
        parser.error("--top-n must be >= 1")
    if args.top_n is not None and not args.sort_field:
        parser.error("--sort-field is required with --top-n")

    source = Path(args.result_manifest).expanduser().resolve()
    rows = read_manifest(source)
    if args.require_success:
        rows = [row for row in rows if row.get("status") == "success"]

    ids = _id_set(args.ids_file, args.ids)
    if ids is not None:
        wanted = set(ids)
        selected = [row for row in rows if row["design_id"] in wanted]
        missing = [item for item in ids if item not in {row["design_id"] for row in selected}]
        if missing:
            raise SystemExit(f"requested design_id not found: {', '.join(missing)}")
    else:
        field = args.sort_field

        numeric: list[tuple[dict[str, str], float]] = []
        text: list[tuple[dict[str, str], str]] = []
        missing: list[dict[str, str]] = []
        for row in rows:
            value = row.get(field, "").strip()
            if not value:
                missing.append(row)
                continue
            try:
                numeric.append((row, float(value)))
            except ValueError:
                text.append((row, value))
        ranked = [
            row for row, _value in sorted(numeric, key=lambda item: item[1], reverse=not args.ascending)
        ]
        ranked.extend(
            row for row, _value in sorted(text, key=lambda item: item[1], reverse=not args.ascending)
        )
        # Missing scores are always last, including for descending selection.
        ranked.extend(missing)
        selected = ranked[: args.top_n]

    output = Path(args.output).expanduser().resolve()
    for row in selected:
        row["source_manifest"] = str(source)
        row["selection_timestamp"] = datetime.now(timezone.utc).isoformat()
    write_tsv(output, selected)
    selection = {
        "source_manifest": str(source),
        "output_manifest": str(output),
        "input_rows": len(rows),
        "output_rows": len(selected),
        "ids": [row["design_id"] for row in selected],
        "require_success": bool(args.require_success),
        "created_at": datetime.now(timezone.utc).isoformat(),
    }
    output.with_name("selection.json").write_text(json.dumps(selection, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"selected {len(selected)}/{len(rows)} rows -> {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
