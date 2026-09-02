function results = run_ml_BC_matched_label_bank()
%RUN_ML_BC_MATCHED_LABEL_BANK Label all still-unopened states in the fixed plane.

results=run_ml_BC_subset_labeling(brusselator_ml_matched_label_bank_config());
end
