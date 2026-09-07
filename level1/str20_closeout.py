"""Compact completed AF3 runs without losing metrics or resume verification."""
import csv
import hashlib
import json
from pathlib import Path

from af3_multichain_archive import atomic_bytes, atomic_json, digest, load_bundle


def checked_path(root, relative):
    path = root / relative
    if path.is_symlink() or root.resolve() not in path.resolve().parents:
        raise ValueError(f'Unsafe closeout path: {path}')
    return path


def verify(root, record):
    if not record['tasks']:
        raise ValueError('Empty closeout manifest')
    for item in record['keep']:
        path = checked_path(root, item['path'])
        # Metadata is moved after the manifest commit; accept either location during recovery.
        if not path.exists() and item.get('before_move'):
            path = checked_path(root, item['before_move'])
        if not path.is_file() or digest(path) != item['sha256']:
            raise ValueError(f'Closeout verification failed: {path}')
    return record


def plan(root):
    tasks = json.loads((root/'tasks.json').read_text())
    statuses = json.loads((root/'status.json').read_text())
    ids = [t['task_id'] for t in tasks]
    if not ids or len(set(ids)) != len(ids) or set(ids) != set(statuses):
        raise ValueError('Task/status manifest mismatch')
    if any(statuses[name]['status'] != 'success' for name in ids):
        raise ValueError('Closeout requires every task to succeed; failed diagnostics retained')
    with (root/'metrics/all.csv').open() as handle:
        rows = list(csv.DictReader(handle))
    if len(rows) != len(ids) or {r['task_id'] for r in rows} != set(ids):
        raise ValueError('Metrics do not cover every task exactly once')
    by_id = {r['task_id']: r for r in rows}
    if (root/'raw').exists() and any(p.is_file() or p.is_symlink() for p in (root/'raw').rglob('*')):
        raise ValueError('Raw outputs remain; finish task-level cleanup first')
    record = dict(schema='str20_compact_v1', state='planned', tasks=tasks, keep=[], delete=[], move=[],
                  input_sha256={}, results={}, merged_logs=[])

    def item(path):
        checked_path(root, path.relative_to(root))
        return dict(path=str(path.relative_to(root)), sha256=digest(path))

    for task in tasks:
        name = task['task_id']
        bundle = checked_path(root, f'archive/{name}')
        archived = load_bundle(bundle)
        if archived['cleanup_state'] != 'complete':
            raise ValueError(f'Unfinished raw cleanup: {name}')
        for key, value in archived['result'].items():
            expected = '' if value is None else str(value)
            if by_id[name].get(key) != expected:
                raise ValueError(f'Metric differs from archive: {name}/{key}')
        record['results'][name] = archived['result']
        record['input_sha256'][name] = digest(bundle/'input.json')
        current_input = root/'inputs'/f'{name}.json'
        if not current_input.is_file() or json.loads(current_input.read_text()) != json.loads((bundle/'input.json').read_text()):
            raise ValueError(f'Input differs from archived prediction: {name}')
        retained = {f'{name}.cif', f'{name}_confidences.json'}
        for path in bundle.iterdir():
            if not path.is_file() or path.is_symlink():
                raise ValueError(f'Unexpected archive entry: {path}')
            record['keep' if path.name in retained else 'delete'].append(item(path))
        if not all((bundle/f).is_file() for f in retained):
            raise ValueError(f'Missing retained structure/confidence: {name}')

    log_bytes = bytearray()
    for path in sorted(root.rglob('*.log')):
        if 'log' in path.relative_to(root).parts:
            continue
        source = item(path)
        data = path.read_bytes()
        if hashlib.sha256(data).hexdigest() != source['sha256']:
            raise ValueError(f'Log changed while merging: {path}')
        log_bytes.extend(f'\n===== {source["path"]} sha256={source["sha256"]} =====\n'.encode())
        source['offset'], source['length'] = len(log_bytes), len(data)
        log_bytes.extend(data)
        log_bytes.extend(b'\n')
        record['merged_logs'].append(source)
        record['delete'].append(dict(path=source['path'], sha256=source['sha256']))
    atomic_bytes(root/'log/run.log', bytes(log_bytes))
    record['keep'].append(item(root/'log/run.log'))
    deleted = {x['path'] for x in record['delete']}
    for path in sorted(root.rglob('*')):
        relative = path.relative_to(root)
        if path.is_symlink():
            raise ValueError(f'Symlink not allowed in closeout: {path}')
        if not path.is_file() or relative.parts[0] in {'archive', 'log'} or str(relative) in deleted:
            continue
        source = item(path)
        if relative.parts[0] in {'inputs', 'jax_cache'}:
            record['delete'].append(source)
        else:
            target = 'log/'+str(relative)
            # Keep old worker lock separate from the stable lock used by the new runner.
            if str(relative) == 'worker.lock':
                target = 'log/legacy_worker.lock'
            record['move'].append(dict(source, target=target))
            record['keep'].append(dict(path=target, sha256=source['sha256'], before_move=source['path']))
    return record


def closeout(root):
    checkpoint = root/'log/closeout.json'
    if checkpoint.exists():
        record = verify(root, json.loads(checkpoint.read_text()))
        if record['state'] == 'complete':
            return record
    else:
        record = plan(root)
        # A recoverable manifest and merged log must be durable before removing any input.
        atomic_json(checkpoint, record)
        verify(root, json.loads(checkpoint.read_text()))
    for item in record['move']:
        source, target = checked_path(root,item['path']), checked_path(root,item['target'])
        if source.exists():
            if digest(source) != item['sha256']:
                raise ValueError(f'Source changed before move: {source}')
            if target.exists():
                raise ValueError(f'Move destination already exists: {target}')
            target.parent.mkdir(parents=True, exist_ok=True)
            source.rename(target)
    verify(root, record)
    merged = (root/'log/run.log').read_bytes()
    for item in record['merged_logs']:
        data = merged[item['offset']:item['offset']+item['length']]
        if hashlib.sha256(data).hexdigest() != item['sha256']:
            raise ValueError(f'Merged log is incomplete: {item["path"]}')
    for item in record['delete']:
        path = checked_path(root,item['path'])
        if path.exists():
            if digest(path) != item['sha256']:
                raise ValueError(f'File changed before cleanup: {path}')
            path.unlink()
    for path in sorted(root.rglob('*'), key=lambda p: len(p.parts), reverse=True):
        if path.is_dir() and not path.is_symlink() and path not in (root/'archive',root/'log'):
            try:
                path.rmdir()
            except OSError:
                pass
    verify(root, record)
    record['state'] = 'complete'
    atomic_json(checkpoint, record)
    return record
