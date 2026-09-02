function status = run_ml_BC_matched_benchmark_pipeline()
%RUN_ML_BC_MATCHED_BENCHMARK_PIPELINE Resume and complete the frozen study.
%   Every scientific stage retains its own refusal-to-overwrite and provenance
%   gates. This wrapper only skips stages whose versioned final artifact exists.

missing_labels=fullfile('experiment_outputs', ...
    'ml_BC_matched_label_bank_missing_v1','ml_BC_matched_missing_labels.mat');
complete_bank=fullfile('experiment_outputs','ml_BC_matched_label_bank_v1', ...
    'complete_label_bank_manifest.json');
frozen_proposals=fullfile('experiment_outputs', ...
    'ml_BC_matched_acquisition_proposals_v1','frozen_proposal_manifest.json');
endpoint_evaluation=fullfile('experiment_outputs', ...
    'ml_BC_matched_acquisition_evaluation_v1', ...
    'matched_acquisition_evaluation_summary.json');
unique_edge_pairs=fullfile('experiment_outputs', ...
    'ml_BC_matched_acquisition_evaluation_v1', ...
    'unique_edge_recovery_manifest.json');
edge_recovery=fullfile('experiment_outputs','ml_BC_matched_edge_recovery_v1', ...
    'ml_BC_edge_recovery.mat');
dynamical_summary=fullfile('experiment_outputs', ...
    'ml_BC_matched_dynamical_benchmark_v1', ...
    'matched_dynamical_benchmark_summary.json');

if ~isfile(missing_labels)
    run_ml_BC_matched_label_bank();
end
if ~isfile(complete_bank)
    run_python('assemble_ml_BC_matched_label_bank.py');
end
if ~isfile(frozen_proposals)
    run_python('freeze_ml_BC_matched_acquisition_proposals.py');
end
if ~isfile(endpoint_evaluation)
    run_python('evaluate_ml_BC_matched_acquisition.py');
end
if ~isfile(unique_edge_pairs)
    run_python('prepare_ml_BC_matched_unique_edge_pairs.py');
end
if ~isfile(edge_recovery)
    run_ml_BC_matched_edge_recovery();
end
if ~isfile(dynamical_summary)
    run_python('analyze_ml_BC_matched_edge_recovery.py');
end

status=struct('missing_labels_complete',isfile(missing_labels), ...
    'complete_bank_assembled',isfile(complete_bank), ...
    'proposals_frozen',isfile(frozen_proposals), ...
    'endpoint_evaluation_complete',isfile(endpoint_evaluation), ...
    'unique_edge_pairs_prepared',isfile(unique_edge_pairs), ...
    'edge_recovery_complete',isfile(edge_recovery), ...
    'dynamical_summary_complete',isfile(dynamical_summary));
disp(struct2table(status,'AsArray',true));
end

function run_python(script)
command=sprintf(['MPLCONFIGDIR=/private/tmp/brusselator_matplotlib ', ...
    'PYTHONPYCACHEPREFIX=/private/tmp/brusselator_pycache ', ...
    'python3 -B "%s"'],script);
[code,output]=system(command,'-echo');
if code~=0
    error('run_ml_BC_matched_benchmark_pipeline:PythonStage', ...
        'Python stage %s failed with status %d:\n%s',script,code,output);
end
end
