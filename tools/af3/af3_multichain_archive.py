"""Validated task bundles and manifest-only cleanup for the AF3 monitor."""
import hashlib
import json
import os
from pathlib import Path

from Bio.PDB import MMCIFParser
from af3_multichain_core import confidence_metrics, protein_ids


def digest(path):
    value = hashlib.sha256()
    with open(path, 'rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            value.update(block)
    return value.hexdigest()


def atomic_bytes(path, content):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(f'.{path.name}.{os.getpid()}.tmp')
    with open(tmp, 'wb') as handle:
        handle.write(content)
        handle.flush()
        os.fsync(handle.fileno())
    tmp.replace(path)


def atomic_json(path, value):
    atomic_bytes(path, (json.dumps(value, allow_nan=False, indent=2) + '\n').encode())


def load_bundle(folder):
    folder = Path(folder)
    with open(folder / 'result.json') as handle:
        record = json.load(handle)
    for item in record['artifacts']:
        path = folder / item['name']
        if path.parent != folder or not path.is_file() or digest(path) != item['sha256']:
            raise ValueError(f'archive verification failed: {path}')
    return record


def process_task(folder_path, chain_ids, cutoff, archive_dir, keep_source, input_dir, noseqs, json_archive_dir=None):
    folder = Path(folder_path).resolve()
    name = folder.name
    try:
        if (folder / 'result.json').exists():
            record = load_bundle(folder)
            current_input = Path(input_dir) / f'{name}.json' if input_dir else None
            if current_input is not None and current_input.is_file():
                saved_input = folder / 'input.json'
                if not saved_input.is_file() or json.loads(saved_input.read_text()) != json.loads(current_input.read_text()):
                    raise ValueError('archive input is missing or differs from current input')
            return record['result']
        summary = next(folder.glob('*_summary_confidences.json'), None)
        structure = next(folder.glob('*model.cif'), None) or next(folder.glob('*.cif'), None)
        if summary is None or structure is None:
            raise ValueError('missing summary or structure')
        detail = summary.with_name(summary.name.replace('_summary_confidences', '_confidences'))
        with open(summary) as handle:
            s = json.load(handle)
        with open(detail) as handle:
            d = json.load(handle)
        parsed = MMCIFParser(QUIET=True).get_structure(name, str(structure))
        if not list(parsed.get_atoms()):
            raise ValueError('empty CIF')
        template = None
        if input_dir and (Path(input_dir) / f'{name}.json').is_file():
            with open(Path(input_dir) / f'{name}.json') as handle:
                template = json.load(handle)
        proteins = protein_ids(template) if template else [
            c.id for c in parsed[0] if any(r.resname in {
                'ALA','ARG','ASN','ASP','CYS','GLN','GLU','GLY','HIS','ILE',
                'LEU','LYS','MET','PHE','PRO','SER','THR','TRP','TYR','VAL'
            } for r in c)
        ]
        result = confidence_metrics(s, d, chain_ids, cutoff, proteins)
        result.update(model_name=name, status='success', error='', protein_chain_order=','.join(proteins))
        if template and not noseqs:
            for entry in template['sequences']:
                if 'protein' in entry:
                    p = entry['protein']
                    for cid in p['id'] if isinstance(p['id'], list) else [p['id']]:
                        result[f'seq_{cid}'] = p['sequence']
        bundle = Path(archive_dir).resolve() / name
        if bundle == folder or folder in bundle.parents:
            raise ValueError('archive must be outside the raw task directory')
        bundle.mkdir(parents=True, exist_ok=True)
        artifacts = []
        sources = [(structure, f'{name}.cif'), (summary, f'{name}_summary_confidences.json'),
                   (detail, f'{name}_confidences.json')]
        if template is not None:
            sources.append((Path(input_dir) / f'{name}.json', 'input.json'))
        for source, filename in sources:
            target = bundle / filename
            atomic_bytes(target, source.read_bytes())
            checksum = digest(source)
            if digest(target) != checksum:
                raise ValueError(f'archive checksum mismatch: {source}')
            artifacts.append({'name': filename, 'sha256': checksum})
            if json_archive_dir and source in (summary, detail):
                atomic_bytes(Path(json_archive_dir) / filename, target.read_bytes())
        result['af3_cif_path'] = str(bundle / f'{name}.cif')
        result['json_archived'] = True
        # Snapshot only files already present; later additions are never deleted.
        cleanup = [{'path': str(p.relative_to(folder)), 'sha256': digest(p)}
                   for p in folder.rglob('*') if p.is_file() and not p.is_symlink()]
        record = {'result': result, 'artifacts': artifacts, 'source': str(folder),
                  'cleanup': cleanup, 'cleanup_state': 'pending'}
        atomic_json(bundle / 'result.json', record)
        load_bundle(bundle)
        return result
    except Exception as exc:
        return {'model_name': name, 'status': 'failed', 'error': str(exc)}


def cleanup_bundle(bundle, raw_root):
    record = load_bundle(bundle)
    source = Path(record['source'])
    if source.parent.resolve() != Path(raw_root).resolve() or source.is_symlink():
        raise ValueError('cleanup source is not a direct raw task directory')
    for item in record['cleanup']:
        path = source / item['path']
        if source.resolve() not in path.resolve().parents or path.is_symlink():
            raise ValueError('unsafe cleanup manifest path')
        if path.is_file() and digest(path) == item['sha256']:
            path.unlink()
    if source.exists():
        for path in sorted(source.rglob('*'), key=lambda p: len(p.parts), reverse=True):
            if path.is_dir() and not path.is_symlink():
                try:
                    path.rmdir()
                except OSError:
                    pass
        try:
            source.rmdir()
        except OSError:
            pass
    record['cleanup_state'] = 'retained_changed_files' if source.exists() else 'complete'
    atomic_json(Path(bundle) / 'result.json', record)
