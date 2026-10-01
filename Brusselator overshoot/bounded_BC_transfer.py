#!/usr/bin/env python3
"""Freeze and evaluate additional transfer tests without changing paper inputs."""
from pathlib import Path
import argparse,hashlib,json,time
from datetime import datetime,timezone
import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist
from probit_gp_kernels import optimize_hyperparameters
from run_ml_BC_acquisition_robustness_proposals import candidate_pairs,common_pair_geometry,greedy_batch,joint_opposite_probability
from freeze_ml_BC_transfer_proposals import rank_method
ROOT=Path(__file__).resolve().parent
OUT=ROOT/'experiment_outputs/bounded_BC_transfer_20261001_v1'
FAMILIES=['flow_cosine_b10','parameter_b10p02']
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def write(p,x):Path(p).write_text(json.dumps(x,indent=2,allow_nan=False)+'\n')
def freeze():
 if OUT.exists():raise FileExistsError(OUT)
 OUT.mkdir()
 cfg=json.loads((ROOT/'ml_bc_transfer_prospective_config.json').read_text())
 protocol=dict(created_utc=datetime.now(timezone.utc).isoformat(),families=FAMILIES,
  family_definitions={'flow_cosine_b10':'Existing physical B/C anchors advanced by autonomous flow for t=1.0; analytic transverse cosine direction, not the old POD3 ribbon.',
    'parameter_b10p02':'b=10.02; Newton/Floquet validation of B,C,E; eight guarded pilot subdivisions of the B-to-C phase-fixed orbit-state chord, no replacement on failure.'},
  flow_advance=1.0,pilot_steps=8,replicates=3,initial_labels=15,endpoint_budget=8,short_pair_max=.2500001,
  guard_steps=10,candidate_duration=220,edge_radius=.02,minimum_shadow=50,period_tolerance=.01,
  design_seeds=[202610011,202610012,202610013],acquisition=cfg,
  design_rule='17x17 grid, anchors 139/151 plus randomized maximin: uniformly choose among 12 farthest unused states, until 15 total; fixed seeds.',
  direct_rule='Shortest normalized opposite-label pair among the same 15 labels; guarded midpoint refinement until separation <=0.25, at most 8 additional labels; stop on any non-B/C outcome.',
  acquisition_rule='Four disjoint pairs per method, 8 endpoint labels; same frozen scores as original prospective experiment; all initial designs and proposal batches frozen before corresponding labels.',
  recovery_rule='Highest-ranked observed B/C proposal only; direct uses its first short pair; ten guarded midpoint steps then independent 220-time trajectory; recovery requires distance<.02, shadow>=50, relative period error<.01.',
  gates='Both B and C, with no third/unresolved class, required in each initial design. Failed families/designs retained in the report. Changed-parameter library names only validated B/C; all nonmatches stay unresolved. No fallback parameter or tuned family.',
  cost_accounting='Report setup/pilot costs separately; per-method initial/acquisition/guard/diagnostic costs charged even when cross-method outcomes are physically shared; measured integration times are descriptive under concurrent jobs.',
  provenance={'base_git_revision':'f084c4ff8670642f63e153623d44a30b573e903c','candidate_labels_precomputed':False,'results_in_manuscript':False})
 deps=['run_bounded_BC_transfer.m','bounded_BC_transfer.py','ml_bc_transfer_prospective_config.json','probit_gp_kernels.py','run_ml_BC_acquisition_robustness_proposals.py','freeze_ml_BC_transfer_proposals.py','brusselator_evolve_and_classify_state.m','brusselator_classify_five_way_trajectory.m','solve_brusselator_1d_forced.m','run_periodic_edge_orbit_newton.m','run_periodic_edge_orbit_floquet.m']
 protocol['implementation_sha256']={p:sha(ROOT/p) for p in deps};write(OUT/'protocol.json',protocol)
 print('Protocol frozen',sha(OUT/'protocol.json'))
def check():
 p=json.loads((OUT/'protocol.json').read_text())
 expected=dict(p['implementation_sha256'])
 amendment=OUT/'execution_amendment.json'
 if amendment.exists():expected.update(json.loads(amendment.read_text())['implementation_sha256'])
 for name,h in expected.items():
  if sha(ROOT/name)!=h:raise RuntimeError('Implementation changed after freeze: '+name)
 return p
def designs(family):
 p=check();fd=OUT/family;G=pd.read_csv(fd/'geometry.csv');X=G[['normalized_alpha_1','normalized_alpha_2']].to_numpy();rows=[]
 if (fd/'designs.csv').exists():raise FileExistsError('Designs already frozen')
 for rep,seed in enumerate(p['design_seeds'],1):
  rng=np.random.default_rng(seed);chosen=[138,150]
  while len(chosen)<15:
   d=cdist(X,X[chosen]).min(axis=1);d[chosen]=-np.inf
   eligible=np.argsort(-d,kind='stable')[:12];chosen.append(int(rng.choice(eligible)))
  rows.extend(dict(replicate=rep,state_index=int(G.iloc[i].state_index),order=j+1) for j,i in enumerate(chosen))
 D=pd.DataFrame(rows);D.to_csv(fd/'designs.csv',index=False)
 pd.DataFrame({'state_index':sorted(D.state_index.unique())}).to_csv(fd/'initial_selection.csv',index=False)
 write(fd/'design_freeze.json',dict(created_utc=datetime.now(timezone.utc).isoformat(),geometry_sha256=sha(fd/'geometry.mat'),csv_sha256=sha(fd/'geometry.csv'),design_sha256=sha(fd/'designs.csv'),protocol_sha256=sha(OUT/'protocol.json')))
 print(family,'initial union',len(D.state_index.unique()))
def proposals(family):
 p=check();fd=OUT/family;freeze=json.loads((fd/'design_freeze.json').read_text());cfg=p['acquisition'];G=pd.read_csv(fd/'geometry.csv');D=pd.read_csv(fd/'designs.csv');ends=set()
 assert sha(fd/'geometry.mat')==freeze['geometry_sha256'];assert sha(fd/'designs.csv')==freeze['design_sha256']
 for rep in range(1,4):
  dest=fd/f'proposals_{rep}.csv'
  if dest.exists():raise FileExistsError(dest)
  started=time.perf_counter();ids=D[D.replicate==rep].state_index.tolist()
  labels=[json.loads((fd/'labels'/f'state_{i:03d}.json').read_text())['final_label'] for i in ids]
  if set(labels)!={'B','C'}:
   write(fd/f'design_{rep}_failed.json',{'status':'initial_binary_gate_failed','labels':labels,'state_indices':ids});continue
  train=G.set_index('state_index').loc[ids];features=cfg['model']['feature_columns'];X=train[features].to_numpy(float);y=np.where(np.array(labels)=='C',1.,-1.);m=cfg['model']
  state,opt=optimize_hyperparameters(X,y,m['initial_starts'],kernel_family=m['kernel_family'],lengthscale_bounds=tuple(m['lengthscale_bounds']),signal_std_bounds=tuple(m['signal_std_bounds']),jitter=m['jitter'],maximum_laplace_iterations=m['maximum_laplace_iterations'],laplace_tolerance=m['laplace_tolerance'],maximum_optimizer_iterations=m['maximum_optimizer_iterations'],function_tolerance=m['function_tolerance'],log_prior_mean=np.log(m['log_prior_mean_parameters']),log_prior_std=np.array(m['log_prior_std']))
  if not state.converged:raise RuntimeError('GP did not converge')
  candidate=G[~G.state_index.isin(ids)].reset_index(drop=True);pairs=candidate_pairs(candidate,cfg);common=common_pair_geometry(state,candidate,pairs,X,cfg)
  q,negative,minimum=joint_opposite_probability(state,common['X'],common['first'],common['second'],cfg['joint_posterior_sampling'],20261001+rep)
  batches=[]
  for method in cfg['methods']:
   batch=greedy_batch(rank_method(method,pairs,common,q,cfg),cfg);batch.insert(0,'method',method);batches.append(batch)
  P=pd.concat(batches,ignore_index=True);P.to_csv(dest,index=False)
  for method,group in P.groupby('method'):
   endpoints=group[['first_state_index','second_state_index']].to_numpy().ravel();assert len(set(endpoints))==8;assert not set(endpoints)&set(ids)
   ends.update(map(int,endpoints))
  write(fd/f'proposal_freeze_{rep}.json',dict(created_utc=datetime.now(timezone.utc).isoformat(),proposal_sha256=sha(dest),design_sha256=sha(fd/'designs.csv'),initial_label_hashes={str(i):sha(fd/'labels'/f'state_{i:03d}.json') for i in ids},GP_and_acquisition_seconds=time.perf_counter()-started,MC_seed=20261001+rep,minimum_joint_eigenvalue=float(minimum),negative_eigenvalue_count=int(negative),candidate_labels_used=False))
 pd.DataFrame({'state_index':sorted(ends)}).to_csv(fd/'endpoints_selection.csv',index=False)
 print(family,'endpoint union',len(ends))
def summarize():
 p=check();rows=[];status={}
 for fam in FAMILIES:
  fd=OUT/fam
  if (fd/'gate_failed.json').exists():status[fam]=json.loads((fd/'gate_failed.json').read_text());continue
  status[fam]=json.loads((fd/'preparation_summary.json').read_text())
  for f in fd.glob('audit_*/summary.json'):rows.append(json.loads(f.read_text()))
  status[fam]['failed_designs']=[json.loads(f.read_text()) for f in fd.glob('design_*_failed.json')]
 frame=pd.DataFrame(rows);frame.to_csv(OUT/'method_results.csv',index=False)
 write(OUT/'family_status.json',status)
 if len(frame):
  summary=frame.groupby(['family','method']).agg(designs=('replicate','count'),short_pairs=('short_pair_found','sum'),recovered=('edge_recovered','sum'),mean_acquisition_labels=('acquisition_labels','mean'),mean_total_labels=('total_label_integrations','mean'),mean_integration_seconds=('total_integration_runtime_seconds','mean'),mean_shadow=('shadow_duration','mean'))
  summary.to_csv(OUT/'summary.csv');print(summary.to_string())
 print(json.dumps(status,indent=2))
if __name__=='__main__':
 a=argparse.ArgumentParser();a.add_argument('action',choices=['freeze','designs','proposals','summary']);a.add_argument('family',nargs='?');args=a.parse_args()
 if args.action=='freeze':freeze()
 elif args.action=='summary':summarize()
 elif args.action=='designs':designs(args.family)
 else:proposals(args.family)
