function results = run_ml_BC_matched_label_bank_shard(shard)
%RUN_ML_BC_MATCHED_LABEL_BANK_SHARD Run one independent frozen label shard.

results=run_ml_BC_subset_labeling( ...
    brusselator_ml_matched_label_bank_shard_config(shard));
end
