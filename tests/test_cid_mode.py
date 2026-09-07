import csv
import json
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from common.manifest import read_manifest
from common.pipeline import prepare_run
from level1.cid import validate_cid_rows
from level1.run_level1 import _context, _run_stage, run_shard
from level1.aggregate_level1 import main as aggregate
from test_cid_multichain import write_task


def make_input(root):
    proteins = [{'protein': {'id': c, 'sequence': 'AAA'}} for c in ('A', 'B')]
    for state in ('apo', 'holo'):
        entities = proteins + ([{'ligand': {'id': 'L', 'ccdCodes': ['ATP']}}] if state == 'holo' else [])
        (root/f'{state}.json').write_text(json.dumps({'dialect': 'alphafold3', 'version': 1, 'sequences': entities}))
    path = root/'input.tsv'
    path.write_text('design_id\tsequence\tchain_ids\tcid_apo_json_path\tcid_holo_json_path\n'
                    'd1\taaa|ggg\tA,B\tapo.json\tholo.json\n')
    return path


class CidModeTests(unittest.TestCase):
    def test_prepare_normalizes_and_resolves_cid_contract(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            run, rows = prepare_run(make_input(root), root/'run', level='level1', shard_size=1, config={'af3_mode': 'cid'})
            self.assertEqual(rows[0]['sequence'], 'AAA:GGG')
            self.assertEqual(rows[0]['cid_apo_json_path'], str(root/'apo.json'))
            self.assertTrue((run/'tasks/shard_00000.tsv').is_file())

    def test_invalid_cid_inputs_fail_before_run_creation(self):
        for field, value in (('chain_ids', 'A'), ('chain_ids', 'A,A'), ('sequence', 'AAA:'),
                             ('sequence', 'AAA:GG*'), ('cid_apo_json_path', ''),
                             ('cid_apo_json_path', 'holo.json')):
            with self.subTest(field=field, value=value), tempfile.TemporaryDirectory() as tmp:
                root = Path(tmp)
                path = make_input(root)
                rows = read_manifest(path)
                rows[0][field] = str(root/value) if value.endswith('.json') else value
                with self.assertRaises(ValueError):
                    validate_cid_rows(rows)

    def test_template_mismatch_and_mixed_systems_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            rows = read_manifest(make_input(root))
            mixed = dict(rows[0], design_id='d2', chain_ids='B,A')
            with self.assertRaisesRegex(ValueError, 'one chain order'):
                validate_cid_rows([rows[0], mixed])
            payload = json.loads((root/'apo.json').read_text())
            payload['sequences'].pop()
            (root/'apo.json').write_text(json.dumps(payload))
            with self.assertRaisesRegex(ValueError, 'protein IDs'):
                prepare_run(root/'input.tsv', root/'run', level='level1', shard_size=1, config={'af3_mode': 'cid'})
            self.assertFalse((root/'run').exists())

    def test_malformed_tsv_rejected(self):
        for contents in ('design_id\tsequence\tsequence\nd\tAAA\tGGG\n',
                         'design_id\tsequence\nd\tAAA\textra\n',
                         'design_id\tsequence\nd\n'):
            with tempfile.TemporaryDirectory() as tmp:
                path = Path(tmp)/'bad.tsv'
                path.write_text(contents)
                with self.assertRaises(ValueError):
                    read_manifest(path)

    def test_cid_default_command_maps_apo_holo_and_chain_order(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            rows = read_manifest(make_input(root))
            context = _context(root, root/'shard', rows, {})
            context['af3_mode'] = 'cid'
            with patch('level1.run_level1.run_logged', return_value=(0, '')) as run:
                _run_stage({}, 'af3', context, root/'log')
            cmd = run.call_args.args[0]
            self.assertEqual(cmd[cmd.index('--template-with-lig')+1], str(root/'holo.json'))
            self.assertEqual(cmd[cmd.index('--template-without-lig')+1], str(root/'apo.json'))
            self.assertEqual(cmd[-3:], ['--sequence-chain-ids', 'A', 'B'])
            self.assertNotIn('--run-esmfold', cmd)

    def test_cid_worker_aggregate_and_resume(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            run, _ = prepare_run(make_input(root), root/'run', level='level1', shard_size=1, config={'af3_mode': 'cid'})
            def fake_stage(config, model, context, log):
                if model == 'esmfold':
                    folder = Path(context['esm_outdir'])
                    folder.mkdir(parents=True)
                    (folder/'d1.pdb').write_text('ATOM      1  CA  ALA A   1       1.000   2.000   3.000  1.00 80.00           C  \nEND\n')
                else:
                    for variant, score in (('with_lig', 0.9), ('without_lig', 0.2)):
                        folder = write_task(Path(context['af3_outdir']), f'd1_af3_{variant}', ['A', 'B'])
                        path = folder/f'd1_af3_{variant}_summary_confidences.json'
                        data = json.loads(path.read_text())
                        data['iptm'] = score
                        path.write_text(json.dumps(data))
                return 0, ''
            with patch('level1.run_level1._run_stage', side_effect=fake_stage) as launch:
                result = run_shard(run, 0)[0]
                self.assertEqual(launch.call_count, 2)
            self.assertEqual(result['af3_holo_iptm'], 0.9)
            self.assertEqual(result['af3_apo_iptm'], 0.2)
            self.assertEqual(result['status'], 'success')
            with patch('level1.run_level1._run_stage') as launch:
                self.assertEqual(run_shard(run, 0)[0]['status'], 'success')
                launch.assert_not_called()
            with patch.object(sys, 'argv', ['aggregate', '--run-dir', str(run)]):
                self.assertEqual(aggregate(), 0)
            for state, score in (('holo', 0.9), ('apo', 0.2)):
                with open(run/f'metrics/cid_{state}.csv') as handle:
                    row = next(csv.DictReader(handle))
                self.assertEqual(float(row['iptm']), score)
                self.assertEqual(row['status'], 'success')
            statuses = (run/'status/status.tsv').read_text()
            self.assertIn('af3_apo', statuses)
            self.assertIn('af3_holo', statuses)


if __name__ == '__main__':
    unittest.main()
