function results = run_physical_ramp_dynamic_edge()
%RUN_PHYSICAL_RAMP_DYNAMIC_EDGE Track the adjacent physical B--A bracket.
cfg = brusselator_physical_ramp_dynamic_edge_config();
results = run_ml_CBA_dynamic_edge_tracking(cfg.dynamic_pair.bracket_id,cfg);
end
