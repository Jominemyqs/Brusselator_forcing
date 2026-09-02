function results = run_ml_BC_transfer_edge_recovery()
%RUN_ML_BC_TRANSFER_EDGE_RECOVERY Guard and track prospective transfer pairs.
results = run_ml_BC_edge_recovery( ...
    brusselator_ml_BC_transfer_edge_recovery_config());
end
