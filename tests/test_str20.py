import csv
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from level1 import str20
from level1 import str20_closeout
from test_cid_multichain import write_task


class Str20Tests(unittest.TestCase):
    def test_three_state_archive_resume_and_failed_diagnostics(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            tasks = []
            for state,ligands in str20.STATES.items():
                name = 's01_'+state
                tasks.append(dict(task_id=name,design_id='s01',state=state))
                payload = dict(name=name,sequences=[{'protein':{'id':c,'sequence':'A'}} for c in 'AB'] +
                    [{'ligand':{'id':c,'ccdCodes':[ccd]}} for c,ccd in zip('CD',ligands)])
                str20.atomic_json(root/'inputs'/f'{name}.json',payload)
            str20.atomic_json(root/'tasks.json',tasks)
            str20.atomic_json(root/'status.json',{t['task_id']:{'status':'pending'} for t in tasks})
            str20.atomic_bytes(root/'jax_cache/test-cache',b'cache')
            str20.atomic_bytes(root/'worker.log',b'worker log contents\n')
            def predict(cmd,**kwargs):
                class Result:
                    returncode = 0
                if '-c' in cmd:
                    return Result()
                name = Path(next(x.split('=',1)[1] for x in cmd if x.startswith('--json_path='))).stem
                state = name.split('_')[1]
                write_task(root/'raw',name,list('ABCD'[:2+len(str20.STATES[state])]))
                return Result()
            with patch.object(str20,'gpu_idle',return_value=True), patch.object(str20.subprocess,'run',side_effect=predict):
                str20.run(root,1)
            with (root/'log/metrics/all.csv').open() as handle:
                self.assertEqual(len(list(csv.DictReader(handle))),3)
            self.assertFalse((root/'raw').exists())
            self.assertFalse((root/'inputs').exists())
            self.assertFalse((root/'jax_cache').exists())
            self.assertEqual({p.name for p in root.iterdir()},{'archive','log'})
            self.assertIn(b'worker log contents', (root/'log/run.log').read_bytes())
            with patch.object(str20,'gpu_idle',return_value=True), patch.object(str20.subprocess,'run',side_effect=predict) as runner:
                str20.run(root,1)
                runner.assert_not_called()  # Compact resume never initializes CUDA.
            for state in str20.STATES:
                self.assertTrue((root/'log/metrics'/f'{state}.csv').exists())
                name='s01_'+state
                self.assertEqual({p.name for p in (root/'archive'/name).iterdir()},
                                 {name+'.cif',name+'_confidences.json'})
            (root/'archive/s01_apo/s01_apo.cif').write_text('corrupted')
            with patch.object(str20.subprocess,'run') as runner:
                with self.assertRaisesRegex(ValueError,'verification failed'):
                    str20.run(root,1)
                runner.assert_not_called()

    def ready(self,root):
        name='s01_apo'
        task=dict(task_id=name,state='apo')
        payload={'sequences':[{'protein':{'id':'A','sequence':'A'}}]}
        str20.atomic_json(root/'inputs'/f'{name}.json',payload)
        str20.atomic_json(root/'tasks.json',[task])
        str20.atomic_json(root/'status.json',{name:dict(status='success')})
        folder=write_task(root/'raw',name,['A'])
        result=str20.process_task(folder,None,15,root/'archive',True,root/'inputs',False)
        str20.cleanup_bundle(root/'archive'/name,root/'raw')
        str20.table(root/'metrics/all.csv',[dict(task,**result)])
        str20.atomic_bytes(root/'jax_cache/test',b'cache')
        str20.atomic_bytes(root/'logs/task.log',b'complete diagnosis\n')

    def test_closeout_manifest_write_failure_preserves_sources(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);self.ready(root)
            with patch.object(str20_closeout,'atomic_json',side_effect=OSError('disk full')):
                with self.assertRaises(OSError): str20_closeout.closeout(root)
            self.assertTrue((root/'inputs/s01_apo.json').exists())
            self.assertTrue((root/'archive/s01_apo/result.json').exists())
            self.assertTrue((root/'logs/task.log').exists())

    def test_interrupted_cleanup_resumes_without_original_inputs(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);self.ready(root)
            original=Path.unlink
            def interrupted(path,*args,**kwargs):
                if path.name=='test': raise OSError('interrupted cache deletion')
                return original(path,*args,**kwargs)
            with patch.object(Path,'unlink',interrupted):
                with self.assertRaises(OSError): str20_closeout.closeout(root)
            self.assertFalse((root/'inputs/s01_apo.json').exists())
            str20.run(root,1)
            self.assertEqual(json.loads((root/'log/closeout.json').read_text())['state'],'complete')
            self.assertFalse((root/'jax_cache').exists())

    def test_failed_or_unverified_run_cannot_be_compacted(self):
        for failure in ['status','metrics','archive']:
            with self.subTest(failure=failure),tempfile.TemporaryDirectory() as tmp:
                root=Path(tmp);self.ready(root)
                if failure=='status': str20.atomic_json(root/'status.json',{'s01_apo':{'status':'failed'}})
                elif failure=='metrics': (root/'metrics/all.csv').write_text('task_id\n')
                else: (root/'archive/s01_apo/s01_apo_confidences.json').write_text('{}')
                with self.assertRaises(ValueError): str20_closeout.closeout(root)
                self.assertTrue((root/'inputs/s01_apo.json').exists())
                self.assertTrue((root/'jax_cache/test').exists())

    def test_metrics_write_failure_preserves_raw(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)
            name='s01_apo'
            payload={'sequences':[{'protein':{'id':'A','sequence':'A'}}]}
            str20.atomic_json(root/'inputs'/f'{name}.json',payload)
            str20.atomic_json(root/'tasks.json',[dict(task_id=name,state='apo')])
            str20.atomic_json(root/'status.json',{name:dict(status='running',prediction_exit_code=0,prediction_input=payload)})
            folder=write_task(root/'raw',name,['A'])
            with patch.object(str20,'gpu_idle',return_value=True), patch.object(str20.subprocess,'run') as runner, patch.object(str20,'table',side_effect=OSError('disk full')):
                runner.return_value.returncode=0
                str20.run(root,1)
            self.assertTrue(folder.exists())
            self.assertEqual(json.loads((root/'status.json').read_text())[name]['status'],'failed')


if __name__=='__main__':
    unittest.main()
