import subprocess,time,json,sys,concurrent.futures
from pathlib import Path
repo=Path(__file__).resolve().parent
fam=sys.argv[1];fd=repo/'experiment_outputs/bounded_BC_transfer_20261001_v1'/fam
matlab='/Applications/MATLAB_R2026a.app/bin/matlab'
def py(action):subprocess.run([str(repo/'.venv/bin/python'),'bounded_BC_transfer.py',action,fam],cwd=repo,check=True)
def run(action,part=1,total=1):
 log=Path('/tmp')/f'bounded_{fam}_{action}_{part}.log'
 with log.open('w') as out:
  subprocess.run([matlab,'-singleCompThread','-batch',f"run_bounded_BC_transfer('{action}','{fam}',{part},{total});"],cwd=repo,stdout=out,stderr=subprocess.STDOUT,check=True)
 print('completed',fam,action,part,flush=True)
while not (fd/'geometry.mat').exists() and not (fd/'gate_failed.json').exists():time.sleep(2)
# geometry.mat is saved before geometry.csv and preparation summary; wait for completion marker.
while not (fd/'preparation_summary.json').exists() and not (fd/'gate_failed.json').exists():time.sleep(2)
if (fd/'gate_failed.json').exists():
 print('Family gate failed:',(fd/'gate_failed.json').read_text(),flush=True);sys.exit(0)
if not (fd/'designs.csv').exists():py('designs')
with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:list(pool.map(lambda part:run('initial',part,2),[1,2]))
if not (fd/'endpoints_selection.csv').exists():py('proposals')
with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:list(pool.map(lambda part:run('endpoints',part,2),[1,2]))
reps=[rep for rep in range(1,4) if (fd/f'proposals_{rep}.csv').exists()]
with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:list(pool.map(lambda rep:run('audit',rep,1),reps))
print('FAMILY_COMPLETE',fam,flush=True)
