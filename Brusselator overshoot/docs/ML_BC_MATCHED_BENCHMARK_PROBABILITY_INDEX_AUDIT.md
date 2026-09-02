# Probability-index audit of the matched B--C acquisition benchmark

## Finding

The frozen primary proposal implementation
`freeze_ml_BC_matched_acquisition_proposals.py` assigns

```python
probability = predict_proba(state, X)[1]
```

but the preserved interface in `probit_gp.py` returns

```text
(probability_positive, latent_mean, latent_variance).
```

Index 1 is therefore the latent mean, not `p(C)`. The implementation hash in
`frozen_implementation_manifest_v7.json` matches the current preserved script,
so this is a property of the frozen run rather than a later source change.

The generated proposal table provides an independent data-level check. Among
144 proposal rows, 96 contain at least one value in the fields
`first_probability_C` or `second_probability_C` outside `[0,1]`; values range
approximately from `-4.90` to `6.10`.

## Scope of impact

Affected:

- marginal endpoint entropies in the pointwise-uncertainty comparator;
- the local probability-gradient orientation factor used by pointwise and
  pair-aware acquisition;
- saved fields described as endpoint probabilities.

Not affected:

- the 289-state PDE label bank;
- the joint latent opposite-sign probability `q`, which is computed from
  `predict_latent_joint`;
- space-filling endpoint novelty;
- endpoint outcome opening;
- guarded PDE refinement, edge tracking, Newton convergence, Floquet analysis,
  or the independently verified periodic orbit `E_BC`.

## Corrective action

The versioned robustness protocol uses index 0 for `p(C)` and derives both
entropy and local-normal orientation from the corrected probability. It leaves
all frozen v1 files unchanged. Its matched replay results supersede the v1
pair-aware-versus-pointwise bracket-yield comparison for manuscript inference.

Until the corrected results are complete, the manuscript should not quote the
v1 `32/48`, `16/48`, and `4/48` comparison as its primary ML evidence. The
dynamical edge-state result remains valid and independent of this acquisition
audit.
