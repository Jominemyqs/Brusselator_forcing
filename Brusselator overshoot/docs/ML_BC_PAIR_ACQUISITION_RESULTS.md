# Joint-posterior pair acquisition around the verified B--C separator

## Outcome

The pair-aware acquisition stage is complete on the retrospective B--C
method-development plane. Four disjoint endpoint pairs were frozen before their
PDE outcomes were generated. Three of the four pairs were genuine B--C
brackets; the remaining pair was C--C.

All three valid brackets completed the unchanged 10-step full-state bisection
without encountering A, R5, R6, ambiguity, or a final-time unresolved
obstruction. Each resulting candidate entered the neighborhood of the
independently Newton-converged periodic edge orbit `E_BC`.

| Frozen method | New endpoint runs | Valid B--C brackets | Exact `E_BC` recoveries |
|---|---:|---:|---:|
| Pointwise entropy acquisition | 8 | 0/4 | 0/4 |
| Budget-matched space filling | 8 | 1/4 | 1/4 |
| Joint-posterior pair acquisition | 8 | 3/4 | 3/4 |

With only four proposals per method, this table is descriptive rather than a
significance test. Its scientific value is that the primary metric is not
classification accuracy: it is recovery of an invariant object that was
computed independently of the ML model.

## Pair-aware exact-edge diagnostics

The exact reference period is `T_E = 4.0029149135788`.

| Pair | Endpoint outcomes | Minimum orbit distance | Edge shadow duration | Direct period | Relative period error |
|---|---|---:|---:|---:|---:|
| `pair_234_250` | C/B | `1.5726e-4` | `150.5` | `3.9969788` | `1.4829e-3` |
| `pair_059_075` | C/B | `7.8706e-5` | `162.5` | `3.9966798` | `1.5576e-3` |
| `pair_180_182` | B/C | `6.6567e-5` | `146.5` | `4.0102821` | `1.8404e-3` |

The neighborhood gate was `2e-2`, so the observed minima are more than two
orders of magnitude below the gate. The candidate final labels remain
`unresolved` at the finite 220-unit horizon; here that is expected and useful,
because the candidates spend 146.5--162.5 time units shadowing an unstable edge
orbit instead of rapidly approaching B or C.

## Acquisition definition

For a candidate pair `(z_i,z_j)`, the primary term is the joint posterior
probability that the latent GP values have opposite signs,

`q_ij = P[f(z_i) f(z_j) < 0 | D]`.

The implemented score multiplies this quantity by a short-distance penalty and
a local-normal orientation factor. Four disjoint pairs are then selected with a
midpoint-diversity penalty. The outcome labels of the proposed endpoints, the
exact edge orbit, and the held-out C--B--A interval near `lambda = 0.276` were
not used in acquisition.

Hyperparameter stress tests used the saved nominal model, a wide-bound MLE, and
a weakly regularized wide-bound MAP fit. All three selected exactly the same
four-pair batch. Nominal versus wide MLE had top-20 Jaccard overlap `1.0`;
nominal versus weak MAP had `0.90476`. This addresses the immediate concern that
the acquisition ranking might be an artifact of a bound-saturated GP fit.

## Cost

The eight endpoint labels used about `862.9 s` of recorded PDE runtime. The 30
bisection classifications used about `2925.6 s`; the complete downstream batch,
including the three 220-unit candidate integrations and diagnostics, used about
`3219.1 s`. Endpoint labeling plus complete downstream recovery therefore used
about `4082.0 s` in this run.

## Scope and next gate

This is strong controlled method-development evidence, not yet a held-out
discovery result. The B--C plane was chosen precisely because its organizing
edge orbit was already known and could serve as a dynamical ground truth.

The next experiment must freeze the pair-selection rule and success criteria,
then expose the untouched narrow C--B--A challenge near `lambda = 0.276` only
once. Success there means more than detecting classification uncertainty: the
method should identify the hidden B window, propose useful adjacent
opposite-outcome brackets, and recover reproducible boundary dynamics without
using the previously known labels to tune the method.

Machine-readable outputs are in:

- `experiment_outputs/ml_BC_pair_acquisition_sensitivity_v1`
- `experiment_outputs/ml_BC_pair_acquisition_frozen_v1`
- `experiment_outputs/ml_BC_pair_acquisition_labeling_v1`
- `experiment_outputs/ml_BC_pair_acquisition_edge_recovery_v1`
- `experiment_outputs/ml_BC_pair_acquisition_benchmark_v1`
