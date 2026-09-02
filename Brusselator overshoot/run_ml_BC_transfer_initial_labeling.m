function results = run_ml_BC_transfer_initial_labeling()
%RUN_ML_BC_TRANSFER_INITIAL_LABELING Label the frozen transfer seed design.
cfg = brusselator_ml_BC_transfer_initial_labeling_config();
manifest = jsondecode(fileread(cfg.ml_subset.prerequisite_file));
if ~manifest.protocol_frozen || manifest.candidate_outcomes_used || ...
        manifest.exact_edge_orbit_used || manifest.initial_state_count~=15
    error('run_ml_BC_transfer_initial_labeling:ProtocolGate', ...
        'The prospective transfer protocol is not safely frozen.');
end
results = run_ml_BC_subset_labeling(cfg);
end
