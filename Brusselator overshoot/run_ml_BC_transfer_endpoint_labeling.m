function results = run_ml_BC_transfer_endpoint_labeling()
%RUN_ML_BC_TRANSFER_ENDPOINT_LABELING Reveal only frozen selected endpoints.
cfg = brusselator_ml_BC_transfer_endpoint_labeling_config();
manifest = jsondecode(fileread(cfg.ml_subset.prerequisite_file));
if ~manifest.proposal_batches_frozen || ...
        manifest.candidate_outcomes_used_for_proposals || ...
        manifest.endpoint_outcomes_used_for_proposals || ...
        manifest.exact_edge_orbit_used
    error('run_ml_BC_transfer_endpoint_labeling:ProposalGate', ...
        'The transfer proposal batch is not safely frozen.');
end
if manifest.unique_selected_endpoint_count~=cfg.ml_subset.expected_unique_count
    error('run_ml_BC_transfer_endpoint_labeling:EndpointCount', ...
        'Frozen endpoint count changed.');
end
results = run_ml_BC_subset_labeling(cfg);
end
