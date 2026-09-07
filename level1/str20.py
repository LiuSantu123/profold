"""Reproducible STR CID three-state AF3 validation on one local GPU."""
import argparse
import csv
import fcntl
import io
import json
import math
import os
import socket
from pathlib import Path
import subprocess
import sys
import time
import traceback
from contextlib import ExitStack

SHARED = Path(os.environ.get('STR20_SOURCE_DIR', '/xcfhome/yhliu/14_magpcr/for_wj/3_filter_v2'))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'tools' / 'af3'))
from af3_multichain_core import fill_template
from af3_multichain_archive import atomic_bytes, atomic_json, process_task, cleanup_bundle, load_bundle
try:
    from .str20_closeout import closeout
except ImportError:
    from str20_closeout import closeout

AF3_PYTHON = os.environ.get('AF3_PYTHON', '/xcfhome/yhliu/002_software/001_conda/envs/af3/bin/python')
AF3_SCRIPT = os.environ.get('AF3_SCRIPT', '/xcfhome/yhliu/002_software/055_alphafold3/AF3_latest_20260904/run_alphafold.py')
AF3_MODEL_DIR = os.environ.get('AF3_MODEL_DIR', '/xcfhome/pubdata/folding/alphafold3/models')
AF3_DB_DIR = os.environ.get('AF3_DB_DIR', '/xcfhome/pubdata/folding/alphafold3')
STATES = {'apo': [], 'holo': ['STR'], 'quat': ['STR', 'A8W']}


def table(path, rows, delimiter=','):
    fields = list(dict.fromkeys(k for row in rows for k in row))
    buf = io.StringIO()
    writer = csv.DictWriter(buf, fields, delimiter=delimiter)
    writer.writeheader()
    writer.writerows(rows)
    atomic_bytes(path, buf.getvalue().encode())


def prepare(root):
    if root.exists():
        raise ValueError(f'Run directory already exists: {root}; use run to resume')
    shortlist = SHARED / 'cid_experiment_shortlist_20260417/STR/STR_selected_summary.csv'
    pool = SHARED / 'done/str_final_v2_abc/str_final_v2_cid_results.csv'
    selected, sequences, families = [], set(), set()

    def add(row, source, short):
        name = row['design_name']
        family = name.rsplit('_', 2)[0]
        a, b = (row['sequence_A'], row['sequence_B']) if short else (row['with_lig_seq_A'], row['with_lig_seq_B'])
        if (a, b) in sequences or (not short and family in families):
            return
        if not a or not b or set(a+b) - set('ACDEFGHIKLMNPQRSTVWY'):
            raise ValueError(f'Invalid sequence: {name}')
        if not short and (a != row['without_lig_seq_A'] or b != row['without_lig_seq_B']):
            raise ValueError(f'Historical state sequence mismatch: {name}')
        apo = float(row['apo_ab_iptm'] if short else row['without_lig_A-B_iptm'])
        holo = float(row['holo_ab_iptm'] if short else row['with_lig_A-B_iptm'])
        if not all(math.isfinite(v) for v in (apo, holo)):
            raise ValueError(f'Invalid historical scores: {name}')
        selected.append(dict(design_id=f's{len(selected)+1:02d}', original_id=name,
                             sequence_A=a, sequence_B=b, family=family, source_csv=str(source),
                             selection='experimental_shortlist' if short else 'delta_iptm_rank_unique_family',
                             historical_apo_iptm=apo, historical_holo_iptm=holo,
                             historical_delta_iptm=holo-apo))
        families.add(family)
        sequences.add((a,b))

    for row in csv.DictReader(shortlist.open()):
        add(row, shortlist, True)
    rows = list(csv.DictReader(pool.open()))
    def score(row):
        try:
            apo, holo = float(row['without_lig_A-B_iptm']), float(row['with_lig_A-B_iptm'])
            return (holo-apo, holo) if math.isfinite(apo+holo) else (-math.inf, -math.inf)
        except ValueError:
            return (-math.inf, -math.inf)
    for row in sorted(rows, key=lambda r: (-score(r)[0], -score(r)[1], r['design_name'])):
        if len(selected) == 20:
            break
        if score(row)[0] > -math.inf:
            add(row, pool, False)
    assert len(selected) == 20
    template = json.loads((SHARED/'template/cid_withlig_str.json').read_text())
    tasks, manifest = [], []
    for state, ligands in STATES.items():
        t = dict(template)
        t['sequences'] = [x for x in template['sequences'] if 'protein' in x]
        t['sequences'] += [{'ligand': {'id': cid, 'ccdCodes': [ccd]}} for cid,ccd in zip('CD',ligands)]
        atomic_json(root/'templates'/f'{state}.json', t)
    for row in selected:
        did = row['design_id']
        manifest.append(dict(design_id=did, sequence=row['sequence_A']+':'+row['sequence_B'], chain_ids='A,B',
                             cid_apo_json_path='templates/apo.json', cid_holo_json_path='templates/holo.json'))
        for state, ligands in STATES.items():
            name = f'{did}_{state}'
            t = json.loads((root/'templates'/f'{state}.json').read_text())
            data = fill_template(t,name,[row['sequence_A'],row['sequence_B']],['A','B'],seed=6028)
            atomic_json(root/'inputs'/f'{name}.json',data)
            tasks.append(dict(task_id=name,design_id=did,original_id=row['original_id'],state=state,
                              chains=','.join('ABCD'[:2+len(ligands)])))
    table(root/'selected.tsv',selected,'\t')
    table(root/'input.tsv',manifest,'\t')
    atomic_json(root/'tasks.json',tasks)
    atomic_json(root/'status.json',{x['task_id']:{'status':'pending'} for x in tasks})
    atomic_json(root/'settings.json',dict(seed=6028,num_diffusion_samples=5,num_recycles=10,
        msa='none',initial_guess=False,python=AF3_PYTHON,script=AF3_SCRIPT))
    atomic_bytes(root/'README.md', (
        '# STR20 three-state validation\n\n'
        '20 historical STR protein pairs; apo=A+B, holo=A+B+STR, quat=A+B+STR+A8W.\n'
        'CCD identities follow historical templates. A8W is an additional ligand, not a protein.\n'
        'First 10: experiment shortlist in original order. Next 10: descending historical\n'
        'holo-minus-apo AB ipTM, then holo ipTM, then ID; exclude duplicate pairs and existing families.\n'
        'selected.tsv records exact source CSVs, sequences and historical scores. This enriched\n'
        'selection validates execution and state handling; it is not an unbiased accuracy benchmark.\n'
        '60 AF3 inputs, seed 6028, 5 diffusion samples each, 10 recycles, no MSA/templates/initial guess.\n'
        'input.tsv is the standard v37 two-state CID input. tasks.json adds the independent quat state.\n'
        'This runner validates the AF3 branch and shared monitor, not the entire ESMFold+AF3 Level1.\n'
        'During prediction each task has a diagnostic log and verified archive bundle.\n'
        'After all tasks succeed, archive retains only CIF and full confidence per task;\n'
        'inputs and JAX cache are deleted. Metrics and provenance move to log/, task logs\n'
        'merge into log/run.log; log/closeout.json verifies compact archives for resume.\n'
        'Protein-ligand and ligand-ligand ipSAE are undefined and left empty.\n'
        'Resume: run this script with run --root THIS_DIRECTORY --gpu GPU --wait-idle.\n'
        'Failed tasks retain raw outputs and logs; three consecutive failures stop the queue.\n'
    ).encode())
    print(f'Prepared {len(selected)} pairs / {len(tasks)} tasks: {root}',flush=True)


def gpu_idle(gpu):
    query = subprocess.run(['nvidia-smi','-i',str(gpu),
        '--query-compute-apps=pid','--format=csv,noheader,nounits'],capture_output=True,text=True)
    if query.returncode or query.stdout.strip():
        return False
    query = subprocess.run(['nvidia-smi','-i',str(gpu),
        '--query-gpu=memory.used,utilization.gpu','--format=csv,noheader,nounits'],capture_output=True,text=True)
    if query.returncode:
        return False
    try:
        memory, utilization = [int(v.strip()) for v in query.stdout.split(',')]
        return memory < 256 and utilization == 0
    except ValueError:
        return False


def run(root, gpu, wait_idle=False):
    with run_lock(root):
        if (root/'log/closeout.json').exists():
            closeout(root)
            return
        statuses = json.loads((root/'status.json').read_text())
        if not statuses or not all(v['status']=='success' for v in statuses.values()):
            run_locked(root,gpu,wait_idle)
            statuses = json.loads((root/'status.json').read_text())
        if statuses and all(v['status']=='success' for v in statuses.values()):
            closeout(root)


def run_lock(root):
    stack = ExitStack()
    try:
        (root/'log').mkdir(exist_ok=True)
        paths = [root/'log/worker.lock']
        if (root/'worker.lock').exists():
            paths.append(root/'worker.lock')
        for path in paths:
            lock = stack.enter_context(path.open('a'))
            fcntl.flock(lock,fcntl.LOCK_EX | fcntl.LOCK_NB)
        return stack
    except BaseException:
        stack.close()
        raise


def run_locked(root, gpu, wait_idle):
    tasks = json.loads((root/'tasks.json').read_text())
    atomic_bytes(root/'worker.pid',f'{os.getpid()}\n'.encode())
    env = dict(os.environ, CUDA_VISIBLE_DEVICES=str(gpu), OMP_NUM_THREADS='8',
               XLA_PYTHON_CLIENT_PREALLOCATE='false', PYTHONUNBUFFERED='1',
               JAX_COMPILATION_CACHE_DIR=str(root/'jax_cache'))
    statuses = json.loads((root/'status.json').read_text())
    results = {}
    (root/'logs').mkdir(exist_ok=True)
    worker = dict(pid=os.getpid(),host=socket.gethostname(),gpu=gpu)
    if wait_idle:
        consecutive_idle = 0
        while consecutive_idle < 3:
            idle = gpu_idle(gpu)
            consecutive_idle = consecutive_idle+1 if idle else 0
            atomic_json(root/'worker.json',dict(worker,status='waiting_gpu',updated=time.time(),
                                              consecutive_idle=consecutive_idle))
            print(f'{time.ctime()} WAIT GPU {gpu}: idle={idle}, checks={consecutive_idle}/3',flush=True)
            if consecutive_idle < 3:
                time.sleep(30)
    elif not gpu_idle(gpu):
        raise RuntimeError('GPU is busy or unavailable; use --wait-idle')
    with (root/'logs/preflight.log').open('a') as log:
        code = subprocess.run([AF3_PYTHON,'-c',
            'import jax; from alphafold3.constants import chemical_components; '
            'print(jax.devices("gpu")); c=chemical_components.Ccd(); '
            'assert "STR" in c and "A8W" in c'],env=env,stdout=log,stderr=subprocess.STDOUT).returncode
    if code:
        atomic_json(root/'worker.json',dict(worker,status='blocked_cuda',updated=time.time()))
        raise RuntimeError('CUDA/CCD preflight failed; see logs/preflight.log')
    atomic_json(root/'worker.json',dict(worker,status='running',updated=time.time()))

    def publish():
        rows = [results[t['task_id']] for t in tasks if t['task_id'] in results]
        if rows:
            table(root/'metrics/all.csv',rows)
        for state in STATES:
            subset = [r for r in rows if r['state']==state]
            if subset:
                table(root/'metrics'/f'{state}.csv',subset)

    def finish(task):
        name = task['task_id']
        bundle = root/'archive'/name
        folder = bundle if (bundle/'result.json').exists() else root/'raw'/name
        result = process_task(folder,None,15.0,root/'archive',True,root/'inputs',False)
        if result['status'] != 'success':
            raise ValueError(result['error'])
        results[name] = dict(task,**result)
        publish()
        cleanup_bundle(bundle,root/'raw')
        statuses[name] = {'status':'success','cleanup':load_bundle(bundle)['cleanup_state'],'updated':time.time()}
        atomic_json(root/'status.json',statuses)

    consecutive_failures = 0
    for task in tasks:
        name = task['task_id']
        try:
            if (root/'archive'/name/'result.json').exists():
                finish(task)
                continue
            # A successful predictor exit is persisted before archive work, allowing crash recovery.
            current_input = json.loads((root/'inputs'/f'{name}.json').read_text())
            if statuses[name].get('prediction_exit_code') == 0 and statuses[name].get('prediction_input') != current_input:
                raise ValueError('Input changed after successful prediction; use a new run directory')
            if statuses[name].get('prediction_exit_code') != 0:
                previous = root/'raw'/name
                if previous.exists():
                    diagnostics = root/'diagnostics'
                    diagnostics.mkdir(exist_ok=True)
                    previous.rename(diagnostics/f'{name}_{time.time_ns()}')
                statuses[name] = {'status':'running','started':time.time()}
                atomic_json(root/'status.json',statuses)
                (root/'logs').mkdir(exist_ok=True)
                cmd = [AF3_PYTHON,AF3_SCRIPT,f'--json_path={root}/inputs/{name}.json',
                       f'--output_dir={root}/raw',f'--model_dir={AF3_MODEL_DIR}',
                       f'--db_dir={AF3_DB_DIR}','--run_data_pipeline=false',
                       '--flash_attention_implementation=xla','--num_diffusion_samples=5','--num_recycles=10']
                print(f'{time.ctime()} START {name}',flush=True)
                with (root/'logs'/f'{name}.log').open('a') as log:
                    log.write(json.dumps(cmd)+'\n'); log.flush()
                    code = subprocess.run(cmd,env=env,stdout=log,stderr=subprocess.STDOUT).returncode
                statuses[name].update(prediction_exit_code=code,prediction_input=current_input)
                atomic_json(root/'status.json',statuses)
                if code:
                    raise RuntimeError(f'AF3 exit {code}; see logs/{name}.log')
            finish(task)
            consecutive_failures = 0
            print(f'{time.ctime()} SUCCESS {name}',flush=True)
        except Exception as exc:
            statuses[name].update(status='failed',error=str(exc),updated=time.time())
            atomic_json(root/'status.json',statuses)
            traceback.print_exc()
            consecutive_failures += 1
            if consecutive_failures >= 3:
                atomic_json(root/'worker.json',dict(worker,status='stopped_failures',updated=time.time()))
                raise RuntimeError('Stopping after three consecutive failures; diagnostics retained')
    atomic_json(root/'worker.json',dict(worker,status='finished' if all(v['status']=='success' for v in statuses.values()) else 'finished_with_failures',updated=time.time(),
                                      counts={s:sum(v['status']==s for v in statuses.values()) for s in ['success','failed','pending']}))
    print('QUEUE FINISHED',flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__, epilog=(
        'Successful runs close out automatically: retain CIF/full confidence in archive/, '
        'move metrics and metadata to log/, merge logs, and delete inputs/JAX cache. '
        'Resume of a compact run verifies hashes without starting CUDA.'))
    parser.add_argument('action',choices=['prepare','run','closeout'],
                        help='prepare inputs; run/resume predictions; or compact a completed run')
    parser.add_argument('--root',type=Path,required=True,help='Run directory (retained after closeout)')
    parser.add_argument('--gpu',type=int,default=1,help='GPU index for new predictions only (default: 1)')
    parser.add_argument('--wait-idle',action='store_true',help='Require three idle checks before CUDA preflight')
    args = parser.parse_args()
    if args.action=='prepare':
        prepare(args.root.resolve())
    elif args.action=='closeout':
        with run_lock(args.root.resolve()):
            record = closeout(args.root.resolve())
            print(f'Closeout {record["state"]}: {len(record["tasks"])} tasks')
    else:
        run(args.root.resolve(),args.gpu,args.wait_idle)
