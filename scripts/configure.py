#!/usr/bin/env python3
"""Resolve machine paths without evaluating shell code or runtime placeholders."""
import argparse
import json
from pathlib import Path
from string import Template
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from common.runner import atomic_json


def render(value, paths):
    if isinstance(value, str):
        return Template(value).substitute(paths)
    if isinstance(value, list):
        return [render(item, paths) for item in value]
    if isinstance(value, dict):
        return {key: render(item, paths) for key, item in value.items()}
    return value


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--template', type=Path, default=ROOT / 'environments/site.example.json')
    parser.add_argument('--conda-root', type=Path, required=True)
    parser.add_argument('--software-root', type=Path, required=True)
    parser.add_argument('--data-root', type=Path, required=True)
    parser.add_argument('--output', type=Path, default=ROOT / 'config.local.json')
    args = parser.parse_args()
    if args.output.exists():
        parser.error('Output exists; use a new --output path to preserve local configuration')
    paths = {'V37_ROOT': str(ROOT), 'CONDA_ROOT': str(args.conda_root.expanduser().resolve()),
             'SOFTWARE_ROOT': str(args.software_root.expanduser().resolve()),
             'DATA_ROOT': str(args.data_root.expanduser().resolve())}
    config = render(json.loads(args.template.read_text()), paths)
    atomic_json(args.output, config)
    print(args.output.resolve())


if __name__ == '__main__':
    main()
