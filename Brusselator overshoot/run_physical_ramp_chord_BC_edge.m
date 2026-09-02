function results = run_physical_ramp_chord_BC_edge()
%RUN_PHYSICAL_RAMP_CHORD_BC_EDGE Test the local chord B--C mechanism.
cfg = brusselator_physical_ramp_BC_dynamic_edge_config();
results = run_ml_CBA_dynamic_edge_tracking(cfg.dynamic_pair.bracket_id,cfg);
end
