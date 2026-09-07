import json
from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from scripts.configure import render
from common.runner import command_template
from level1.cid import cid_command


class PortabilityTests(unittest.TestCase):
    def test_profile_preserves_runtime_placeholders_and_spaces(self):
        roots = {name: '/tmp/with spaces/' + name for name in
                 ('V37_ROOT', 'CONDA_ROOT', 'SOFTWARE_ROOT', 'DATA_ROOT')}
        config = render(json.loads((ROOT / 'environments/site.example.json').read_text()), roots)
        cmd = command_template(config['commands']['protenix'],
                               {'protenix_input': '/tmp/a b.json', 'protenix_outdir': '/tmp/out'})
        self.assertEqual(cmd[-2], '/tmp/a b.json')
        self.assertEqual(cmd[1], 'PROTENIX_PYTHON=/tmp/with spaces/CONDA_ROOT/envs/protenix/bin/python')
        self.assertNotIn('${', json.dumps(config))

    def test_missing_site_variable_fails(self):
        with self.assertRaises(KeyError):
            render('${MISSING}/bin/python', {})

    def test_cid_uses_bundled_scripts(self):
        cmd = cid_command({})
        self.assertEqual(Path(cmd[1]), ROOT / 'tools/af3/af3_single_node_v18.py')
        self.assertTrue(Path(cmd[cmd.index('--monitor-script') + 1]).is_file())
