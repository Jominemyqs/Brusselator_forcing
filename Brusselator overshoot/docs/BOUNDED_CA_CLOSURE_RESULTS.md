# Bounded closure of the evolved C--A candidate

## Prespecified stopping rule

The closure experiment was deliberately bounded:

1. refine the sampled C/A interval `[0.875,0.9]` with the unchanged five-way
   A/B/C/R5/R6 classifier;
2. stop immediately if any midpoint reaches B, R5, R6, or remains unresolved
   at the maximum adaptive horizon;
3. run dynamic C--A rebracketing only if five clean C/A midpoint labels reduce
   the interval to width `7.8125e-4`.

This rule was fixed before the first new midpoint outcome was computed.

## Result

The first midpoint gives

```text
alpha = 0.8875
outcome = B
tested frozen time = 200
decision time = 116
late median B-orbit distance = 2.5544e-3
```

The state is not marginal under the existing classifier: its late distances to
A and C are approximately `0.5460` and `0.3809`, respectively. The experiment
therefore terminates with status `obstructed_by_B` after one new PDE trajectory.
No dynamic C--A edge tracking is run.

At the currently sampled coordinates, the local chord ordering is now

```text
B through alpha = 0.825
C at alpha = 0.85 and 0.875
B at alpha = 0.8875
A from alpha = 0.9
```

Thus the apparent C--A interval contains another B-basin intrusion. It yields
candidate C--B and B--A intervals `[0.875,0.8875]` and `[0.8875,0.9]`, each of
width `0.0125`, but they are not claimed to be direct adjacencies.

## Scientific interpretation

This result completes the bounded closure and strengthens the multiclass
interpretation:

> Along this one evolved state-space chord, successive resolution exposes the
> sampled ordering B--C--B--A rather than a direct B--A or C--A transition.

One-dimensional endpoint labels therefore cannot safely define binary basin
boundaries in this region. A smooth folded basin tongue can produce this
ordering, so the result does **not** establish fractal, riddled, or Wada basin
geometry. Establishing any of those stronger properties would require a
separate multidirectional and resolution-aware basin study.

The known B--C mechanism remains the only locally verified separator mechanism:
the periodic edge orbit `E_BC` has exact Newton/Floquet and unstable-branch
support. No C--A or B--A edge state is claimed.

## Reproducibility

- Configuration: `brusselator_ml_CA_closure_config.m`
- Driver: `run_ml_CA_guarded_refinement.m`
- Output: `experiment_outputs/ml_heldout_CBA_evolved_CA_refinement_v1`

The output records the raw trajectory, complete physical distance history,
classification decision, source-state provenance, configuration, metadata,
summary tables, and the stopped bracket. Completed output is never overwritten.
