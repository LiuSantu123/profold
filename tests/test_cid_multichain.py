import argparse
import json
import math
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import gemmi
import numpy as np

SCRIPTS = Path(__file__).resolve().parents[1] / 'tools' / 'af3'
sys.path.insert(0, str(SCRIPTS))
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from af3_multichain_core import confidence_metrics, directional_ipsae, fill_template, read_designs
from af3_multichain_archive import cleanup_bundle, load_bundle, process_task
from af3_single_node_v18 import batch_generate_json_cid, run_af3_prediction
from af3_monitor_v16 import check_design_completion, monitor_loop
from level1.run_level1 import _write_af3_fasta
from common.model_parsers import parse_af3


def fixture(chains):
    n = len(chains)
    return ({'iptm': 0.8, 'ptm': 0.9, 'chain_ptm': [0.1 * (i + 1) for i in range(n)],
             'chain_iptm': [0.2] * n,
             'chain_pair_iptm': [[(i * n + j) / (n * n) for j in range(n)] for i in range(n)],
             'chain_pair_pae_min': [[1] * n for _ in range(n)]},
            {'token_chain_ids': chains, 'atom_chain_ids': chains, 'atom_plddts': [80] * n,
             'pae': [[float(i * n + j + 1) for j in range(n)] for i in range(n)]})


def write_task(root, name, chains):
    folder = root / name
    folder.mkdir(parents=True)
    summary, detail = fixture(chains)
    (folder / f'{name}_summary_confidences.json').write_text(json.dumps(summary))
    (folder / f'{name}_confidences.json').write_text(json.dumps(detail))
    pdb = ''.join(f'ATOM  {i:5d}  CA  ALA {c}   1       1.000   2.000   3.000  1.00 80.00           C  \n'
                  for i, c in enumerate(chains, 1)) + 'END\n'
    structure = gemmi.read_pdb_string(pdb)
    structure.make_mmcif_document().write_file(str(folder / f'{name}_model.cif'))
    return folder


class MultichainTests(unittest.TestCase):
    def test_wrapper_resume_checks_archived_input(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            inputs = root/'inputs'
            inputs.mkdir()
            name = 'd_af3_with_lig'
            payload = {'name': name, 'sequences': [{'protein': {'id': 'A', 'sequence': 'A'}}]}
            source = inputs/f'{name}.json'
            source.write_text(json.dumps(payload))
            folder = write_task(root/'raw', name, ['A'])
            result = process_task(folder, [], 15, root/'archive', False, inputs, False)
            self.assertEqual(result['status'], 'success', result)
            cleanup_bundle(root/'archive'/name, root/'raw')
            args = argparse.Namespace(archive_dir=str(root/'archive'), temp_dir=tmp)
            with patch('af3_single_node_v18._run_af3_prediction') as runner:
                self.assertEqual(run_af3_prediction(args, inputs, root/'raw'), 0)
                runner.assert_not_called()
                payload['sequences'][0]['protein']['sequence'] = 'GG'
                source.write_text(json.dumps(payload))
                with self.assertRaisesRegex(ValueError, 'differs'):
                    run_af3_prediction(args, inputs, root/'raw')
                runner.assert_not_called()
            stale = process_task(root/'archive'/name, [], 15, root/'archive', False, inputs, False)
            self.assertEqual(stale['status'], 'failed')

    def test_failed_task_preserves_diagnostics(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            folder = write_task(root/'raw', 'd_af3_with_lig', ['A', 'B'])
            (folder/'d_af3_with_lig_confidences.json').unlink()
            (folder/'error.log').write_text('diagnostic')
            result = process_task(folder, [], 15, root/'archive', False, None, True)
            self.assertEqual(result['status'], 'failed')
            self.assertTrue((folder/'error.log').is_file())
            self.assertTrue((folder/'d_af3_with_lig_model.cif').is_file())

    def test_keep_source_and_csv_failure_prevent_cleanup(self):
        for keep_source in (False, True):
            with self.subTest(keep_source=keep_source), tempfile.TemporaryDirectory() as tmp:
                root = Path(tmp)
                raw = root/'raw'
                for variant in ('with_lig', 'without_lig'):
                    write_task(raw, f'd_af3_{variant}', ['A'])
                args = argparse.Namespace(noseqs='false', run_esm_compare='false', esm_results_json=None,
                    archive_dir=str(root/'archive'), af3_output_dir=str(raw), csv_output=str(root/'metrics.csv'),
                    processes=1, chain_ids=[], pae_cutoff=15, keep_source=keep_source, json_input_dir=None,
                    json_archive_dir=None, usalign_path=None, once=True)
                if keep_source:
                    monitor_loop(args)
                else:
                    with patch('af3_monitor_v16.save_results', side_effect=OSError('disk full')):
                        with self.assertRaises(OSError):
                            monitor_loop(args)
                self.assertTrue((raw/'d_af3_with_lig/d_af3_with_lig_model.cif').is_file())

    def test_one_to_four_chains_and_v37_fasta(self):
        for n in range(1, 5):
            with self.subTest(n=n), tempfile.TemporaryDirectory() as tmp:
                path = Path(tmp) / 'in.fa'
                _write_af3_fasta([{'design_id': 'd', 'sequence': '|'.join(['AAA'] * n)}], path, 'for_wj_cid')
                designs = read_designs(path)
                self.assertEqual(len(designs[0]) - 1, n)
                chains = list('ABCD'[:n])
                s, d = fixture(chains)
                metrics = confidence_metrics(s, d, protein_chains=chains)
                self.assertEqual(len([k for k in metrics if k.endswith('_ipae') and '-' in k]), n * (n - 1) // 2)
                if n == 1:
                    self.assertIsNone(metrics['overall_ipae'])

    def test_template_nonstandard_ids_and_no_truncation(self):
        template = {'sequences': [{'protein': {'id': ['X', 'Y'], 'sequence': 'OLD'}},
                                  {'ligand': {'id': 'L', 'ccdCodes': ['ATP']}}]}
        result = fill_template(template, 'd', ['AAA', 'GGG'], ['Y', 'X'])
        self.assertEqual(result['sequences'][0]['protein']['sequence'], 'GGG')
        self.assertEqual(result['sequences'][1]['protein']['sequence'], 'AAA')
        self.assertEqual(result['sequences'][2], template['sequences'][1])
        self.assertEqual(result['sequences'][0]['protein']['unpairedMsa'], '')
        with self.assertRaises(ValueError):
            fill_template(template, 'd', ['AAA'])
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / 'a.json').write_text(json.dumps(template))
            (root / 'b.json').write_text(json.dumps({'sequences': template['sequences'][:1]}))
            batch_generate_json_cid([('d', 'AAA', 'GGG')], root/'a.json', root/'b.json', root/'out')
            a = json.loads((root/'out/d_af3_with_lig.json').read_text())
            b = json.loads((root/'out/d_af3_without_lig.json').read_text())
            self.assertEqual(a['modelSeeds'], b['modelSeeds'])

    def test_chain_subset_indices_and_asymmetry(self):
        s, d = fixture(['X', 'Y', 'Z'])
        result = confidence_metrics(s, d, ['Y', 'Z'], protein_chains=['X', 'Y', 'Z'])
        self.assertAlmostEqual(result['Y_ptm'], 0.2)
        self.assertAlmostEqual(result['Y-Z_iptm'], (5 / 9 + 7 / 9) / 2)
        self.assertEqual(result['Y-Z_ipae'], 7.0)
        self.assertEqual(result['Y-Z_pae_min'], 1.0)
        self.assertNotEqual(result['Y-Z_ipsae_forward'], result['Y-Z_ipsae_reverse'])
        with self.assertRaises(ValueError):
            confidence_metrics(s, d, ['missing'])
        d['pae'] = [[1]]
        with self.assertRaises(ValueError):
            confidence_metrics(s, d)

    def test_d0res_formula_and_ligand_na(self):
        block = np.ones((2, 40)) * 5
        expected_d0 = max(1.0, 1.24 * (40 - 15) ** (1 / 3) - 1.8)
        self.assertAlmostEqual(directional_ipsae(block, 15), 1 / (1 + (5 / expected_d0) ** 2))
        self.assertEqual(directional_ipsae(block, 1), 0)
        s, d = fixture(['A', 'L'])
        self.assertIsNone(confidence_metrics(s, d, protein_chains=['A'])['A-L_ipsae'])

    def test_completion_without_esm_and_failed_variant(self):
        tasks = {'d': {v: {'status': 'success', 'protein_chain_order': 'A,B,C,D'}
                       for v in ('with_lig', 'without_lig')}}
        self.assertEqual(check_design_completion('d', {}, tasks), (True, '2/2', []))
        self.assertFalse(check_design_completion('d', {}, tasks, True)[0])
        tasks['d']['without_lig']['status'] = 'failed'
        self.assertFalse(check_design_completion('d', {}, tasks)[0])

    def test_archive_cleanup_resume_and_v37_metrics(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            raw, archive = root/'raw', root/'raw/archive'
            for variant in ('with_lig', 'without_lig'):
                write_task(raw, f'd_af3_{variant}', ['A', 'B', 'C', 'D'])
            args = argparse.Namespace(noseqs='false', run_esm_compare='false', esm_results_json=None,
                archive_dir=str(archive), af3_output_dir=str(raw), csv_output=str(root/'metrics.csv'),
                processes=1, chain_ids=[], pae_cutoff=15, keep_source=False, json_input_dir=None,
                json_archive_dir=None, usalign_path=None, once=True)
            monitor_loop(args)
            self.assertFalse((raw/'d_af3_with_lig').exists())
            self.assertIn('with_lig_C-D_ipae', (root/'metrics.csv').read_text())
            metrics = parse_af3(raw, 'd', mode='for_wj_cid')
            self.assertEqual(metrics['af3_status'], 'success')
            self.assertIn('af3_with_lig_C-D_ipsae', metrics)
            self.assertEqual(metrics['af3_holo_C-D_ipsae'], metrics['af3_with_lig_C-D_ipsae'])
            self.assertIn('af3_apo_C-D_ipae', metrics)
            (root/'metrics.csv').unlink()
            monitor_loop(args)
            self.assertIn('with_lig_C-D_ipae', (root/'metrics.csv').read_text())
            (archive/'d_af3_with_lig/d_af3_with_lig.cif').write_text('corrupt')
            self.assertNotEqual(parse_af3(raw, 'd', mode='for_wj_cid')['af3_status'], 'success')

    def test_bad_archive_prevents_cleanup_and_new_files_survive(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            raw, archive = root/'raw', root/'archive'
            folder = write_task(raw, 'd_af3_with_lig', ['A'])
            result = process_task(folder, [], 15, archive, False, None, True)
            self.assertEqual(result['status'], 'success', result)
            (folder/'new.log').write_text('new diagnostic')
            cleanup_bundle(archive/folder.name, raw)
            self.assertTrue((folder/'new.log').exists())
            (archive/folder.name/f'{folder.name}.cif').write_text('corrupt')
            with self.assertRaises(ValueError):
                cleanup_bundle(archive/folder.name, raw)
            self.assertTrue((folder/'new.log').exists())


if __name__ == '__main__':
    unittest.main()
