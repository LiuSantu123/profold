#!/usr/bin/env python3
"""Read-only dependency/path checks; never starts inference or downloads models."""
import argparse
import importlib.util
import json
from pathlib import Path
import shutil
import sys

ROOT = Path(__file__).resolve().parents[1]


def check(config, level):
    failures = []

    def record(label, ok):
        print(f'{"OK" if ok else "MISSING"} {label}')
        if not ok:
            failures.append(label)

    for name in ('numpy', 'Bio', 'gemmi', 'pandas'):
        record(f'controller module {name}', importlib.util.find_spec(name) is not None)
    for name in ('af3_multichain_core.py', 'af3_multichain_archive.py', 'af3_monitor_v16.py', 'af3_single_node_v18.py'):
        record(f'bundled {name}', (ROOT / 'tools/af3' / name).is_file())

    def path(value, executable=False):
        p = Path(value).expanduser()
        ok = bool(shutil.which(str(p))) if executable else p.exists()
        record(str(p), ok)

    commands = config.get('commands', {})
    models = {'level1': ('esmfold', 'af3'), 'level2': ('protenix', 'boltz2', 'opendde'),
              'level3': ('netsolp', 'temberture')}[level]
    keys = {'esmfold': ('esmfold_bin',), 'af3': ('af3_python', 'af3_script', 'af3_model_dir', 'af3_db_dir'),
            'boltz2': ('boltz_bin',), 'protenix': ('protenix_wrapper',),
            'opendde': ('opendde_wrapper',), 'netsolp': ('netsolp_wrapper',),
            'temberture': ('temberture_python', 'temberture_script')}
    for model in models:
        cmd = commands.get(model, commands.get('boltz') if model == 'boltz2' else None)
        if cmd:
            import shlex
            tokens = shlex.split(cmd) if isinstance(cmd, str) else cmd
            path(tokens[0], executable=True)
            for token in tokens[1:]:
                value = token.split('=', 1)[-1]
                if value.startswith('/') and '{' not in value:
                    path(value)
            print(f'NOTE {model}: custom command paths checked; CLI/model compatibility needs a GPU smoke test')
        else:
            for key in keys[model]:
                record(f'config {key}', bool(config.get(key)))
                if config.get(key):
                    path(config[key], executable=key.endswith(('_python', '_bin')))
    if level == 'level1' and config.get('af3_mode') == 'cid':
        for key, default in [('cid_wrapper', 'af3_single_node_v18.py'), ('cid_monitor', 'af3_monitor_v16.py')]:
            path(config.get(key, str(ROOT / 'tools/af3' / default)))
    print('Path checks only: CUDA, weights and model versions are not validated.')
    return 1 if failures else 0


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config', type=Path, required=True)
    parser.add_argument('--level', choices=['level1', 'level2', 'level3'], required=True)
    args = parser.parse_args()
    sys.exit(check(json.loads(args.config.read_text()), args.level))
