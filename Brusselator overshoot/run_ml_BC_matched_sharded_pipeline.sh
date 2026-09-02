#!/bin/zsh
set -euo pipefail

matlab_bin=/Applications/MATLAB_R2026a.app/bin/matlab
log_dir=experiment_outputs/ml_BC_matched_label_bank_shards_v1
for shard in 1 2 3 4; do
  log_file=${log_dir}/shard_0${shard}_unsandboxed_matlab.log
  "${matlab_bin}" -batch "run_ml_BC_matched_label_bank_shard(${shard})" > "${log_file}" 2>&1
done

MPLCONFIGDIR=/private/tmp/brusselator_matplotlib PYTHONPYCACHEPREFIX=/private/tmp/brusselator_pycache python3 -B assemble_ml_BC_matched_label_bank.py
MPLCONFIGDIR=/private/tmp/brusselator_matplotlib PYTHONPYCACHEPREFIX=/private/tmp/brusselator_pycache python3 -B freeze_ml_BC_matched_acquisition_proposals.py
MPLCONFIGDIR=/private/tmp/brusselator_matplotlib PYTHONPYCACHEPREFIX=/private/tmp/brusselator_pycache python3 -B evaluate_ml_BC_matched_acquisition.py
PYTHONPYCACHEPREFIX=/private/tmp/brusselator_pycache python3 -B prepare_ml_BC_matched_unique_edge_pairs.py
"${matlab_bin}" -batch "run_ml_BC_matched_edge_recovery"
PYTHONPYCACHEPREFIX=/private/tmp/brusselator_pycache python3 -B analyze_ml_BC_matched_edge_recovery.py
