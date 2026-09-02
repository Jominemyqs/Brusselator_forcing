function cfg = brusselator_ml_BC_transfer_endpoint_labeling_config()
%BRUSSELATOR_ML_BC_TRANSFER_ENDPOINT_LABELING_CONFIG Selected endpoint labels.

cfg = brusselator_ml_BC_transfer_initial_labeling_config();
cfg.experiment_name = 'ml_BC_transfer_endpoint_labeling_v1';
proposal_dir = fullfile('experiment_outputs','ml_BC_transfer_proposals_v1');
endpoint_file = fullfile(proposal_dir,'frozen_transfer_endpoint_manifest.csv');
if ~isfile(endpoint_file)
    error('brusselator_ml_BC_transfer_endpoint_labeling_config:MissingProposals', ...
        'Transfer proposals must be frozen before endpoint labeling.');
end
selected = readtable(endpoint_file);
cfg.ml_subset.kind = 'prospective_transfer_selected_endpoints';
cfg.ml_subset.state_index_file = endpoint_file;
cfg.ml_subset.expected_unique_count = height(selected);
cfg.ml_subset.required_design_roles = {'acquisition_pool'};
cfg.ml_subset.prerequisite_file = fullfile(proposal_dir, ...
    'frozen_transfer_proposal_manifest.json');
cfg.ml_subset.final_filename = 'ml_BC_transfer_endpoint_labels.mat';
cfg.ml_subset.final_csv = 'endpoint_labels.csv';
cfg.ml_subset.summary_csv = 'endpoint_labeling_summary.csv';
cfg.ml_subset.figure_filename = 'endpoint_labels.png';
cfg.ml_subset.interpretation = [ ...
    'actual PDE labels only for the union of endpoints selected by the ', ...
    'three hash-frozen transfer acquisition batches'];
end
