function cfg = brusselator_physical_ramp_BC_dynamic_edge_config()
%BRUSSELATOR_PHYSICAL_RAMP_BC_DYNAMIC_EDGE_CONFIG Track local chord B--C.
%   The B/C pair was discovered by multiclass-safe rebracketing of the
%   adjacent physical B/A forcing-return states. This test asks whether the
%   pair converges to the independently known periodic edge orbit E_BC.

cfg = brusselator_physical_ramp_dynamic_edge_config();
cfg.experiment_name = 'physical_ramp_chord_BC_dynamic_edge_v1';
cfg.dynamic_pair.bracket_id = 'physical_ramp_chord_BC';
cfg.dynamic_pair.left_source_file = fullfile('experiment_outputs', ...
    'physical_ramp_BA_dynamic_edge_v1','raw_cases','event01_step08.mat');
cfg.dynamic_pair.right_source_file = fullfile('experiment_outputs', ...
    'physical_ramp_BA_dynamic_edge_v1','raw_cases','event01_step10.mat');
cfg.dynamic_pair.left_coordinate = 0.26171875;
cfg.dynamic_pair.right_coordinate = 0.2626953125;
cfg.dynamic_pair = rmfield(cfg.dynamic_pair, ...
    {'left_coordinate_column','right_coordinate_column'});
cfg.dynamic_pair.left_outcome = 'B';
cfg.dynamic_pair.right_outcome = 'C';
cfg.dynamic_pair.exact_edge_comparison = true;
cfg.dynamic_pair.maximum_rebisection_steps = 10;
cfg.dynamic_pair.segment_maximum_duration = 220;
cfg.dynamic_pair.claim_scope = [ ...
    'one-event B--C dynamic edge tracking on the frozen chord locally ', ...
    'connecting adjacent physical B/A forcing-return states at N=400 and ', ...
    'Lx=40; recovery of E_BC establishes nearby frozen separator geometry, ', ...
    'not direct forcing access to C or proof that E_BC causes the B/A ', ...
    'protocol transition'];
end
