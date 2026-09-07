"""Shared AF3 chain mapping and confidence metrics (token-aligned)."""
import copy
import itertools
import math
import re

import numpy as np


def read_designs(path):
    records, name, chunks = [], None, []
    def append():
        if name is None:
            return
        chains = ''.join(chunks).upper().replace('|', ':').split(':')
        if not all(chains):
            raise ValueError(f'{name}: empty protein chain')
        records.append((name, *chains))
    with open(path) as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                append()
                name = line[1:].split()[0]
                if not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]*', name):
                    raise ValueError(f'unsafe design ID: {name}')
                chunks = []
            else:
                if name is None:
                    raise ValueError('sequence before FASTA header')
                chunks.append(''.join(line.split()))
    append()
    if not records or len({r[0] for r in records}) != len(records):
        raise ValueError('empty FASTA or duplicate design IDs')
    return records


def protein_ids(template):
    ids = []
    for entry in template['sequences']:
        if 'protein' in entry:
            value = entry['protein']['id']
            ids.extend(value if isinstance(value, list) else [value])
    if not ids or len(ids) != len(set(ids)):
        raise ValueError('template requires unique protein chain IDs')
    return ids


def fill_template(template, name, sequences, chain_ids=None, seed=1):
    ids = protein_ids(template)
    order = list(chain_ids or ids)
    if len(sequences) != len(order) or len(order) != len(set(order)) or set(order) != set(ids):
        raise ValueError(f'{name}: FASTA chains and template protein IDs do not match: {ids}')
    mapping = dict(zip(order, sequences))
    output = copy.deepcopy(template)
    output['name'], output['modelSeeds'] = name, [seed]
    expanded = []
    for entry in output['sequences']:
        if 'protein' not in entry:
            expanded.append(entry)
            continue
        value = entry['protein']['id']
        for cid in value if isinstance(value, list) else [value]:
            protein = copy.deepcopy(entry['protein'])
            protein.update(id=cid, sequence=mapping[cid], unpairedMsa='', pairedMsa='', templates=[])
            protein.pop('templateMsa', None)
            protein.pop('unpairedMsaPath', None)
            protein.pop('pairedMsaPath', None)
            expanded.append({'protein': protein})
    output['sequences'] = expanded
    return output


def number(value):
    return float(value) if isinstance(value, (int, float)) and math.isfinite(value) else None


def mean(values):
    valid = [v for v in values if number(v) is not None]
    return float(np.mean(valid)) if valid else None


def directional_ipsae(block, cutoff):
    """Dunbrack ipSAE d0res: row-specific valid-residue d0, then best row.

    Formula reference: github.com/DunbrackLab/IPSAE (v4, 2026-01-03).
    Protein-protein pairs only; small-molecule atom tokens are not residues.
    """
    valid = block < cutoff
    counts = valid.sum(axis=1)
    d0 = np.maximum(1.0, 1.24 * (np.maximum(26, counts) - 15) ** (1 / 3) - 1.8)
    scores = np.where(valid, 1 / (1 + (block / d0[:, None]) ** 2), 0).sum(axis=1)
    return float(np.max(scores / np.maximum(counts, 1))) if len(scores) else 0.0


def confidence_metrics(summary, detail, chains=None, cutoff=15.0, protein_chains=None):
    token_ids = detail.get('token_chain_ids', [])
    atom_ids = detail.get('atom_chain_ids', [])
    order = list(dict.fromkeys(token_ids))
    if not order:
        raise ValueError('token_chain_ids is required for chain-resolved metrics')
    target = list(chains or order)
    if len(target) != len(set(target)) or set(target) - set(order):
        raise ValueError('analysis chain IDs must be unique and present in confidence JSON')
    pae = np.asarray(detail.get('pae', []), dtype=float)
    plddt = np.asarray(detail.get('atom_plddts', []), dtype=float)
    if pae.shape != (len(token_ids), len(token_ids)) or not np.isfinite(pae).all() or (pae < 0).any():
        raise ValueError('PAE shape/values do not match token_chain_ids')
    if len(plddt) != len(atom_ids) or not len(plddt) or not np.isfinite(plddt).all():
        raise ValueError('atom confidence values do not match atom_chain_ids')
    out = {key: summary.get(key) for key in ('iptm', 'ptm', 'ranking_score', 'fraction_disordered', 'has_clash')}
    out.update(chain_order=','.join(target), af3_mean_plddt=float(plddt.mean()),
               overall_mean_pae=float(pae.mean()), metrics_schema='multichain_v2',
               ipae_method='bidirectional_mean_cross_token_pae', ipsae_method='d0res_protein_pairs')
    indices = {c: np.flatnonzero(np.array(token_ids) == c) for c in target}
    for c in target:
        pos = order.index(c)
        ai = np.flatnonzero(np.array(atom_ids) == c)
        ti = indices[c]
        out.update({f'{c}_len_tokens': len(ti), f'{c}_len_atoms': len(ai),
                    f'{c}_mean_plddt': float(plddt[ai].mean()) if len(ai) else None,
                    f'{c}_mean_intra_pae': float(pae[np.ix_(ti, ti)].mean())})
        for key in ('ptm', 'iptm'):
            values = summary.get(f'chain_{key}') or []
            out[f'{c}_{key}'] = number(values[pos]) if pos < len(values) else None
    pair_ipaes, pair_scores = [], []
    for a, b in itertools.combinations(target, 2):
        key = f'{a}-{b}'
        i, j = order.index(a), order.index(b)
        for source, label in (('chain_pair_iptm', 'iptm'), ('chain_pair_pae_min', 'pae_min')):
            matrix = summary.get(source) or []
            if matrix and (len(matrix) != len(order) or any(len(row) != len(order) for row in matrix)):
                raise ValueError(f'{source} does not match chain order')
            ab, ba = (number(matrix[i][j]), number(matrix[j][i])) if matrix else (None, None)
            out.update({f'{key}_{label}_forward': ab, f'{key}_{label}_reverse': ba,
                        f'{key}_{label}': mean([ab, ba])})
        ab, ba = pae[np.ix_(indices[a], indices[b])], pae[np.ix_(indices[b], indices[a])]
        out[f'{key}_ipae_forward'], out[f'{key}_ipae_reverse'] = float(ab.mean()), float(ba.mean())
        out[f'{key}_ipae'] = mean([float(ab.mean()), float(ba.mean())])
        pair_ipaes.append(out[f'{key}_ipae'])
        protein_pair = protein_chains is not None and {a, b} <= set(protein_chains)
        forward = directional_ipsae(ab, cutoff) if protein_pair else None
        reverse = directional_ipsae(ba, cutoff) if protein_pair else None
        out.update({f'{key}_ipsae_forward': forward, f'{key}_ipsae_reverse': reverse,
                    f'{key}_ipsae': max(forward, reverse) if protein_pair else None,
                    f'{key}_ipsae_max': max(forward, reverse) if protein_pair else None,
                    f'{key}_ipsae_min': min(forward, reverse) if protein_pair else None,
                    f'{key}_ipsae_mean': mean([forward, reverse])})
        pair_scores.append(out[f'{key}_ipsae'])
    out['overall_ipae'] = mean(pair_ipaes)
    out['overall_pair_ipsae_mean'] = mean(pair_scores)
    out['overall_iptm'] = out.get('iptm') if len(order) > 1 else None
    return out
