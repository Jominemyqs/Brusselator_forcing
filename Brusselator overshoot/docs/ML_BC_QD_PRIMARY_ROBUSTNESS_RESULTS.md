# qD-primary robustness of task-aware B--C edge-bracket acquisition

## Status and purpose

The versioned masked-replay extension `ml_BC_acquisition_qD_primary_v1`
completed on 2026-08-16. It resolves the final-method mismatch between the
corrected development benchmark, whose reference score was `qDO`, and the
prospective physical-chord/POD3 experiment, which froze the more parsimonious
`qD` score.

The extension changes only the declared pair-aware score. It uses the identical:

- 289-state B/C label bank;
- 40 nested matched initial designs;
- initial-label counts `n0=15,25,35`;
- GP model, optimizer, kernels, and joint-posterior samples;
- pair admissibility and four-pair/eight-label budget;
- pointwise and space-filling comparators; and
- paired statistical procedures.

The geometry, complete-bank, and matched-design hashes are byte-for-byte equal
to the corrected `qDO`-primary study. No new PDE integrations were performed.
The old outputs remain immutable and provide the orientation-factor ablation.

## Primary qD result

The reference condition is `n0=25`, squared-exponential kernel, `qD` score, and
distance scale `0.18`.

| Method | Valid brackets | Pooled yield | Designs with at least one | Brackets per endpoint label |
|---|---:|---:|---:|---:|
| Pair-aware `qD` | 111/160 | 0.69375 | 38/40 | 0.346875 |
| Pointwise uncertainty | 56/160 | 0.35000 | 28/40 | 0.175000 |
| Space filling | 0/160 | 0.00000 | 0/40 | 0.000000 |

At the matched-design level, `qD` exceeds pointwise uncertainty by a mean
`1.375` valid brackets among four, with paired-bootstrap 95% interval
`[1.000,1.750]`. The win/tie/loss counts are `31/6/3`; the exact two-sided sign
test after removing ties is `7.66e-7`.

Relative to space filling, the mean advantage is `2.775`, with interval
`[2.400,3.150]` and win/tie/loss counts `38/2/0`.

### Budget-efficiency curve

The fraction of matched designs with at least one valid bracket is:

| Revealed endpoint labels | Pair-aware `qD` | Pointwise | Space filling |
|---:|---:|---:|---:|
| 2 | 0.675 | 0.225 | 0.000 |
| 4 | 0.800 | 0.500 | 0.000 |
| 6 | 0.900 | 0.625 | 0.000 |
| 8 | 0.950 | 0.700 | 0.000 |

Failures remain in the denominator; these are unconditional design-level
discovery probabilities on the fixed replay plane.

## qD sensitivity across benchmark choices

| Scenario | Pair-aware `qD` | Pointwise | Space filling | qD design coverage |
|---|---:|---:|---:|---:|
| `n0=15` | 87/160 | 35/160 | 6/160 | 35/40 |
| Reference `n0=25` | 111/160 | 56/160 | 0/160 | 38/40 |
| `n0=35` | 122/160 | 69/160 | 7/160 | 39/40 |
| Matérn `3/2` | 105/160 | 49/160 | 0/160 | 38/40 |
| Matérn `5/2` | 103/160 | 53/160 | 0/160 | 37/40 |
| Distance scale `0.75` | 94/160 | 56/160 | 0/160 | 35/40 |
| Distance scale `1.25` | 114/160 | 56/160 | 0/160 | 38/40 |

For all seven prespecified `qD` conditions, the paired-bootstrap interval for
`qD` minus pointwise valid-bracket count excludes zero. Mean advantages range
from `0.950` to `1.450` brackets per design. Pair-aware `qD` also exceeds space
filling in every condition.

This establishes internal robustness on one prelocalized PC1--PC2 state family;
it does not establish transfer across physical systems. The separate
physical-chord/POD3 experiment supplies one prospective local transfer test.

## Acquisition ablation and final rule

| Pair score | Valid brackets | Design coverage | Median normalized pair length |
|---|---:|---:|---:|
| `q` | 126/160 | 39/40 | 0.250000 |
| `qD` | 111/160 | 38/40 | 0.176777 |
| `qDO` | 112/160 | 38/40 | 0.176777 |

`qD` and `qDO` are indistinguishable at the replicate level. For `qD-qDO`,
the mean difference is `-0.025`, the bootstrap interval is `[-0.150,0.100]`,
and win/tie/loss counts are `3/33/4`.

Unpenalized `q` has higher raw yield, but it selects the maximum admissible
normalized separation `0.25` for nearly every proposal. The `D` term therefore
implements a deliberate yield--length tradeoff required by the downstream edge
tracking task. The orientation factor `O` adds no detectable bracket-yield
benefit beyond `qD`.

The frozen final-method sequence is consequently:

```text
joint pair probability q establishes task alignment
    -> D regularizes toward shorter edge brackets
    -> O is rejected as unnecessary in the corrected ablation
    -> qD is frozen before prospective transfer
```

The manuscript should use `qD` as the primary method, retain `q` and `qDO` as
ablations, and replace the former `112/160` primary number with `111/160`.

## Relationship to the prospective experiment

The two positive experiments now have complementary roles:

1. **Repeated development benchmark:** `qD` gives `111/160` brackets versus
   `56/160` and `0/160` across 40 masked matched designs on the development
   plane, with positive paired advantages under every tested perturbation.
2. **Frozen prospective transfer:** the already-frozen `qD` rule gives `2/4`
   brackets versus `0/4` and `0/4` using actual new PDE labels on the
   physical-chord/POD3 family; its prespecified bracket recovers `E_BC`.

The held-out C--B--A failure remains unchanged and establishes the scope
boundary between local pairwise bracket construction and unseen rare-basin
discovery.

## Integrity record

Versioned outputs are stored in:

```text
experiment_outputs/ml_BC_acquisition_qD_primary_protocol_v1
experiment_outputs/ml_BC_acquisition_qD_primary_proposals_v1
experiment_outputs/ml_BC_acquisition_qD_primary_evaluation_v1
```

Principal hashes are:

```text
complete label bank
3f1e9fef80b33923faf45b844446b693ee11f5d9fac8af9554cd8d547eebe0c5

nested matched initial designs
8e6246a6d35e427f6c4e2f0314f347fcbed66d274a999cd70c41e92fd8f1fd62

frozen proposal table
43a795d98831e27be284be6579341ffc4e5c7edbb2430bb27eda6dd3c77fa3e7

final summary
f663a7c38d122f4bed96851f81485d7926b2070821098632d2286c61a4d249c1
```

All output hashes recorded in the final summary were recomputed and passed.
The new reference and `qDO` ablation batches exactly reproduce the corresponding
`qD` and `qDO` batches in the prior corrected suite.
