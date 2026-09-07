"""CID manifest contract and paired AF3 command construction."""
from __future__ import annotations

import json
import re
from pathlib import Path

CID_FIELDS = ("design_id", "sequence", "chain_ids", "cid_apo_json_path", "cid_holo_json_path")
SCRIPTS = Path(__file__).resolve().parents[1] / 'tools' / 'af3'


def validate_cid_rows(rows):
    """Validate before creating a run or launching any GPU work."""
    batch = None
    templates = {}
    for row in rows:
        name = row['design_id']
        missing = [key for key in CID_FIELDS if not row.get(key, '').strip()]
        if missing:
            raise ValueError(f'{name}: CID requires {", ".join(missing)}')
        sequence = row['sequence'].strip().upper().replace('|', ':')
        chains = sequence.split(':')
        ids = [c.strip() for c in row['chain_ids'].split(',')]
        if any(not re.fullmatch(r'[ACDEFGHIKLMNPQRSTVWYX]+', seq) for seq in chains):
            raise ValueError(f'{name}: sequence requires nonempty protein chains (20 amino acids or X), separated by colon')
        if len(ids) != len(chains) or len(set(ids)) != len(ids) or any(not re.fullmatch('[A-Z]+', c) for c in ids):
            raise ValueError(f'{name}: chain_ids must be unique uppercase IDs matching the sequence chain count')
        row['sequence'], row['chain_ids'] = sequence, ','.join(ids)
        pair = tuple(row[key] for key in ('cid_apo_json_path', 'cid_holo_json_path'))
        signature = (tuple(ids), pair)
        if batch is not None and signature != batch:
            raise ValueError('CID run requires one chain order and one apo/holo template pair; split different systems into separate runs')
        batch = signature
        for state, path in zip(('apo', 'holo'), pair):
            if path not in templates:
                with open(path) as handle:
                    template = json.load(handle)
                all_ids, proteins, ligand_count = [], [], 0
                for entry in template.get('sequences', []):
                    if len(entry) != 1:
                        raise ValueError(f'{path}: each entity must have one type')
                    kind, entity = next(iter(entry.items()))
                    entity_ids = entity.get('id', [])
                    entity_ids = entity_ids if isinstance(entity_ids, list) else [entity_ids]
                    all_ids.extend(entity_ids)
                    if kind == 'protein':
                        proteins.extend(entity_ids)
                    if kind == 'ligand':
                        ligand_count += len(entity_ids)
                if len(set(all_ids)) != len(all_ids):
                    raise ValueError(f'{path}: duplicate entity chain IDs')
                templates[path] = (set(proteins), ligand_count)
            proteins, ligand_count = templates[path]
            if proteins != set(ids):
                raise ValueError(f'{name}: {state} template protein IDs differ from chain_ids')
            if (state == 'apo' and ligand_count) or (state == 'holo' and not ligand_count):
                raise ValueError(f'{name}: apo must have no ligand entity; holo must have at least one')
    return rows


def cid_command(config):
    command = [str(config.get('af3_python', '/xcfhome/zpzeng/originrepo/alphafold3/alphafold3_env/bin/python')),
               str(config.get('cid_wrapper', SCRIPTS / 'af3_single_node_v18.py'))]
    values = {
        '--fasta': '{af3_fasta}', '--template-with-lig': '{cid_holo_json_path}',
        '--template-without-lig': '{cid_apo_json_path}',
        '--af3-script': config.get('af3_script', '/xcfhome/zpzeng/originrepo/alphafold3/run_alphafold.py'),
        '--af3-env': command[0], '--model-dir': config.get('af3_model_dir', '/xcfhome/pubdata/folding/alphafold3/models'),
        '--db-dir': config.get('af3_db_dir', '/xcfhome/pubdata/folding/alphafold3'),
        '--monitor-script': config.get('cid_monitor', SCRIPTS / 'af3_monitor_v16.py'),
        '--json-output-dir': '{af3_input_dir}', '--af3-output-dir': '{af3_outdir}',
        '--csv-output': '{af3_csv}', '--archive-dir': '{af3_archive_dir}',
        '--json-archive-dir': '{af3_json_archive_dir}', '--temp-dir': '{af3_temp_dir}',
        '--pae-cutoff': config.get('cid_pae_cutoff', 15),
    }
    for flag, value in values.items():
        command.extend([flag, str(value)])
    if config.get('cid_keep_source', False):
        command.append('--keep-source')
    return command
