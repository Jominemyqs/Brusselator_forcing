# Corrected robustness study for B--C edge-bracket acquisition

## Status

The versioned masked-replay study `ml_BC_acquisition_robustness_v1` completed
on 2026-08-15. It used the existing 289-state PDE label bank and generated no
new PDE integrations. The protocol, 40 nested matched designs, nine
one-factor-at-a-time scenarios, GP implementation, proposal code, evaluation
code, and statistical summaries are hash-frozen in:

```text
experiment_outputs/ml_BC_acquisition_robustness_protocol_v1
experiment_outputs/ml_BC_acquisition_robustness_proposals_v1
experiment_outputs/ml_BC_acquisition_robustness_evaluation_v1
```

The complete bank hash is
`3f1e9fef80b33923faf45b844446b693ee11f5d9fac8af9554cd8d547eebe0c5`.

## Question and scope

The study addresses a conditional acquisition problem:

> Given a prelocalized state family containing a candidate pairwise B--C
> separator, how should a fixed endpoint-label budget be allocated to generate
> nearby opposite-outcome pairs that can initialize edge tracking?

It does not establish global basin discovery, rare-basin discovery, transfer to
an independent PDE state family, or realized end-to-end wall-clock savings.
The 40 designs are algorithmic replay replicates on one physical POD plane, not
40 independent PDE experiments.

## Probability-index audit

Before this extension was frozen, the primary matched-v1 proposal script was
found to use the latent mean returned at index 1 of `predict_proba` as if it
were the class probability returned at index 0. In the v1 proposal table, 96 of
144 rows consequently contain a reported endpoint “probability” outside
`[0,1]`.

The mismatch affects the v1 pointwise entropy comparator and the local
probability-gradient orientation factor used in pointwise and pair-aware
scoring. It does not affect:

- the PDE label bank;
- the joint latent opposite-sign probability `q`;
- geometric space filling;
- guarded edge tracking;
- Newton/Floquet computations; or
- the independently verified periodic orbit `E_BC`.

The corrected robustness implementation uses index 0 for `p(C)`. All 4,320
new proposal rows contain probabilities in `[0,1]`. The results below supersede
the v1 `32/48`, `16/48`, `4/48` comparison for manuscript inference. The v1
files remain unchanged as an audit record.

## Frozen design

- State family: existing `17 x 17`, `N=400`, `Lx=40`, frozen-`b=10` B--C
  landing-aware POD plane.
- Label bank: 159 B and 130 C states; no third or unresolved class.
- Matched designs: 40 randomized-maximin ordered designs with fixed B and C
  anchors.
- Nested initial-label prefixes: `n0=15`, `25`, and `35`.
- Kernels: ARD squared exponential, Matérn `3/2`, and Matérn `5/2`.
- Pair-score variants: `q`, `qD`, and `qDO`.
- Soft distance scales: `0.75`, `1.0`, and `1.25` times `0.18`.
- Hard admissible normalized separation: fixed at
  `[0.1249999, 0.2500001]`.
- Batch: four disjoint pairs, eight newly revealed endpoint labels.
- Independent statistical unit: matched initial-design replicate.

Proposal generation required 200 unique GP fits. Fits and joint posterior
samples were reused across scenarios with identical training count and kernel.
The complete proposal phase ran in `44.04 s`; this is algorithmic replay time,
not a PDE cost.

## Primary corrected reference result

The reference condition is `n0=25`, squared-exponential kernel, `qDO` score,
and distance scale `0.18`.

| Method | Valid brackets | Pooled yield | Designs with at least one | Brackets per endpoint label |
|---|---:|---:|---:|---:|
| Pair-aware `qDO` | 112/160 | 0.700 | 38/40 | 0.350 |
| Pointwise uncertainty | 56/160 | 0.350 | 28/40 | 0.175 |
| Space filling | 0/160 | 0.000 | 0/40 | 0.000 |

The pair-aware method therefore doubles reference bracket yield relative to
the corrected pointwise comparator. At the replicate level, its mean advantage
is `1.400` valid brackets among four, with paired-bootstrap 95% interval
`[1.000, 1.800]`. Pair-aware wins, ties, and loses in `31/6/3` designs; the
exact two-sided sign-test value after removing ties is `7.66e-7`.

Relative to space filling, the mean advantage is `2.800` brackets with interval
`[2.400, 3.175]` and win/tie/loss counts `38/2/0`.

### Budget-efficiency curve

The probability of obtaining at least one bracket after each number of revealed
endpoint labels is:

| Revealed labels | Pair-aware | Pointwise | Space filling |
|---:|---:|---:|---:|
| 2 | 0.675 | 0.225 | 0.000 |
| 4 | 0.850 | 0.500 | 0.000 |
| 6 | 0.900 | 0.625 | 0.000 |
| 8 | 0.950 | 0.700 | 0.000 |

These are unconditional design-level discovery probabilities; failures remain
in the denominator.

## Robustness across benchmark perturbations

The table reports valid brackets among 160 proposals, followed in parentheses
by the number of designs with at least one bracket.

| Scenario | Pair-aware `qDO` | Pointwise | Space filling |
|---|---:|---:|---:|
| `n0=15` | 91 (36/40) | 35 (25/40) | 6 (5/40) |
| Reference `n0=25` | 112 (38/40) | 56 (28/40) | 0 (0/40) |
| `n0=35` | 125 (39/40) | 69 (33/40) | 7 (7/40) |
| Matérn `3/2` | 100 (38/40) | 49 (29/40) | 0 (0/40) |
| Matérn `5/2` | 101 (35/40) | 53 (29/40) | 0 (0/40) |
| Distance scale `0.75` | 100 (36/40) | 56 (28/40) | 0 (0/40) |
| Distance scale `1.25` | 114 (38/40) | 56 (28/40) | 0 (0/40) |

For every prespecified `qDO` scenario, the paired-bootstrap interval for
pair-aware minus pointwise valid-bracket count excludes zero. Mean advantages
range from `1.100` to `1.450` brackets per design. Pair-aware also exceeds
space filling under every perturbation.

The absolute space-filling result varies with initial-design size but is exactly
zero under all three `n0=25` kernels and all three distance scales. This should
be described as a property of the fixed replay plane and the declared
short-pair geometry, not as a universal failure of space filling.

## Pair-score ablation

| Pair score | Valid brackets | Pooled yield | Designs with at least one | Brackets per endpoint label |
|---|---:|---:|---:|---:|
| `q` | 126/160 | 0.7875 | 39/40 | 0.39375 |
| `qD` | 111/160 | 0.69375 | 38/40 | 0.346875 |
| `qDO` | 112/160 | 0.7000 | 38/40 | 0.3500 |

`qD` and `qDO` are indistinguishable at the replicate level: their mean
difference is `0.025`, bootstrap interval `[-0.100, 0.150]`, and win/tie/loss
counts `4/33/3`. The orientation factor therefore provides no detectable
valid-bracket gain beyond the joint probability and distance penalty on this
plane.

`q` alone has higher raw bracket yield than either distance-penalized variant.
However, this is not an unconditional dominance result. It chooses normalized
separation `0.25` in 159/160 proposals, while `qD` and `qDO` choose separation
`0.17678` in 141/160 and 140/160 proposals, respectively, with the remainder
mostly at `0.125`. The distance penalty intentionally trades some probability
of straddling the separator for a shorter starting bracket. The appropriate
methodological conclusion is:

> The joint opposite-sign probability supplies the principal task alignment;
> the distance penalty controls a yield--length tradeoff; and the orientation
> heuristic does not improve bracket yield beyond `qD` in this benchmark.

## Dynamical verification and its boundary

This replay extension evaluates endpoint bracket discovery. It does not rerun
guarded PDE refinement for every newly replayed proposal and therefore does not
turn the 40 designs into new dynamical experiments. The independent dynamical
result remains:

- the B--C separator is organized by the Newton-converged periodic orbit
  `E_BC`;
- its single nontrivial unstable Floquet direction connects to B and C;
- the orbit and branch geometry persist across `N=200,400,600`;
- every one of the 20 unique valid brackets selected in the original guarded
  recovery audit reached `E_BC`.

Because the original acquisition orientation field was affected by the index
audit, those 20 recoveries establish the edge-tracking mechanism and show that
valid local proposals can reach it; they are not a conditional dynamical audit
of all newly corrected robustness proposals. Manuscript language should keep
these two evidential roles separate.

The prospective held-out C--B--A result also remains unchanged:

```text
8 frozen labels -> 7 A + 1 C + 0 B.
```

It is the explicit evidence that success in a prelocalized binary family does
not imply discovery of a narrow unseen basin component.

## Frozen manuscript claims

The manuscript may state:

> On a prelocalized B--C state family, corrected pair-aware acquisition doubles
> valid-bracket yield relative to pointwise uncertainty under the reference
> eight-label budget and retains a positive paired advantage across initial
> label counts, three kernel families, and prespecified distance-scale
> perturbations.

It may also state:

> Ablation identifies the joint opposite-sign probability as the principal
> source of task alignment. The short-pair penalty controls a yield--length
> tradeoff, while the orientation factor adds no detectable bracket-yield
> benefit on the controlled plane.

It must not state:

- that the GP or acquisition method discovers arbitrary basin structure;
- that 40 replay designs are independent physical experiments;
- that the study demonstrates realized PDE wall-clock savings;
- that every corrected replay bracket has been dynamically edge-tracked;
- that space filling generally cannot locate separators; or
- that the method can reliably discover a narrow unseen third basin.

## Integrity record

The final summary is

```text
experiment_outputs/ml_BC_acquisition_robustness_evaluation_v1/
    ml_BC_acquisition_robustness_summary.json
```

Its SHA-256 hash is
`49bdeada687040444e7ed669e2272f9812355d4238b404d4be95b751d5fa383e`.
The summary contains hashes for every CSV, both overview figures, the frozen
protocol manifest, the frozen proposal manifest, and the complete PDE label
bank.
