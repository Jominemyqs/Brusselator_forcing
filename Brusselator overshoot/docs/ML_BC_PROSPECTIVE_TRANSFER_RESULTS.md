# Prospective transfer test for task-aware B--C edge-bracket acquisition

## Status

The versioned experiment `ml_BC_transfer_*_v1` completed on 2026-08-15
(2026-08-16 UTC). It is the first prospective end-to-end test of the frozen
pair-aware acquisition rule on a state family distinct from the original
PC1--PC2 development plane.

The complete workflow generated actual frozen-PDE labels, rather than replaying
labels from the 289-state development bank:

```text
freeze geometry and initial design
    -> label 15 initial states with the frozen PDE
    -> freeze three four-pair acquisition batches
    -> label only the selected endpoints with the frozen PDE
    -> open endpoint outcomes
    -> guarded B/C refinement of the prespecified valid proposal
    -> comparison with the independently computed E_BC orbit
```

## Scientific and methodological question

The experiment asks a deliberately conditional transfer question:

> After a candidate B--C region has been prelocalized, can the frozen
> pair-aware rule transfer to a differently parameterized state family and
> produce a dynamically useful edge-tracking bracket more efficiently than
> corrected pointwise uncertainty or space filling under the same labeling
> budget?

It does not ask whether the method can discover an arbitrary basin, an unseen
third class, or an edge region from the full PDE state space without prior
localization.

## Independent state family

The test family is a `17 x 17`, `N=400`, `Lx=40` physical-chord ribbon. Its
longitudinal direction is the chord between independently labeled B and C
states found on the physical ramp-return landing chord. Its transverse
direction is landing-aware POD mode 3 after orthogonalization against that
chord. Thus the test family is not the PC1--PC2 plane used to develop and
replay the acquisition method.

The basis-projector distance from the original PC1--PC2 plane is `1.5152`,
which exceeds the frozen distinctness threshold `0.2`. The source B/C labels
were used only to prelocalize and orient the family. Candidate outcomes, the
development label bank, the held-out C--B--A region, and the exact `E_BC` orbit
were not used to construct the geometry or acquisition batches.

The geometry and a fixed 15-state coarse initial design were hash-frozen before
any nonanchor outcomes were generated. The initial design produced:

| Outcome | Count |
|---|---:|
| B | 11 |
| C | 4 |
| Other known outcome | 0 |
| Unresolved | 0 |

This passed the prespecified binary-family gate.

## Frozen acquisition comparison

All methods used the same 15 labeled states and probabilistic fit. Each method
proposed four disjoint pairs and received the same budget of eight newly
revealed endpoint labels:

1. geometric space filling;
2. corrected pointwise posterior uncertainty; and
3. pair-aware `qD`, combining joint posterior probability of opposite labels
   with a short-pair penalty.

The proposal batches were frozen before endpoint outcomes were generated. The
three batches contained 23 unique physical endpoints because one endpoint was
shared across methods. Their actual PDE classifications were 8 B and 15 C,
with no third or unresolved outcome.

| Method | Valid B/C brackets | Proposals | Designs with a bracket | Endpoint labels to first bracket |
|---|---:|---:|---:|---:|
| Pair-aware `qD` | 2 | 4 | yes | 4 |
| Pointwise uncertainty | 0 | 4 | no | -- |
| Space filling | 0 | 4 | no | -- |

The result is descriptive evidence from one prospective transfer family, not a
new repeated statistical comparison. The 40-design masked-replay robustness
study remains the appropriate source for statistical evidence within the
original prelocalized plane.

## Prespecified dynamical verification

For each method that produced a valid bracket, the frozen rule selected its
highest-ranked valid proposal for guarded refinement. Only pair-aware
acquisition qualified. Its selected proposal was the method's rank-2 pair,
state 165 (B) and state 166 (C), with normalized separation `0.125`.

Ten adaptive five-way midpoint classifications produced only B or C. Tested
horizons increased as the bracket tightened, with the final midpoint requiring
integration to `t=400` before resolving to C. The segment-coordinate width
decreased from `1` to `9.765625e-4` without encountering A, `P5.06`, `P6.58`,
or a terminal unresolved case.

The resulting boundary trajectory passed the independent exact-orbit audit:

| Quantity | Result |
|---|---:|
| Exact-edge neighborhood recovered | yes |
| Minimum phase-aware physical distance to `E_BC` | `1.3142986e-4` |
| Edge-shadow duration | `182.5` |
| Direct period estimate | `3.9972223` |
| Independently computed `E_BC` period | `4.0029149` |
| Relative period error | `1.4221e-3` |

The candidate's final outcome is recorded as unresolved because it remains near
the boundary orbit over the audit window. That field is not an acquisition or
edge-recovery failure; the prespecified success criterion is recovery of the
exact edge neighborhood and its recurrent period.

## Frozen interpretation

The experiment supports the following claim:

> On one prospectively labeled physical-chord/POD3 state family, the frozen
> pair-aware acquisition rule produced two usable B--C brackets under an
> eight-label budget, while corrected pointwise uncertainty and space filling
> produced none. Guarded refinement of the prespecified pair-aware proposal
> recovered the independently established periodic edge orbit `E_BC`.

Together with the repeated masked-replay result, this shows two complementary
forms of evidence: statistical robustness across sparse designs on the
development plane, and one successful prospective end-to-end transfer to a
differently parameterized local family.

The experiment does **not** establish:

- global or full-state basin discovery;
- rare or unseen class discovery;
- statistical generality across independent physical separators;
- discovery of a previously unknown edge orbit, because `E_BC` was retained
  only as an independent post-acquisition audit target; or
- superiority of either comparator in general.

The previously frozen narrow C--B--A failure remains the relevant evidence that
pairwise bracket acquisition does not automatically solve rare-basin discovery.

## Reproducibility record

Versioned outputs are stored in:

```text
experiment_outputs/ml_BC_transfer_geometry_v1
experiment_outputs/ml_BC_transfer_protocol_v1
experiment_outputs/ml_BC_transfer_initial_labeling_v1
experiment_outputs/ml_BC_transfer_proposals_v1
experiment_outputs/ml_BC_transfer_endpoint_labeling_v1
experiment_outputs/ml_BC_transfer_evaluation_v1
experiment_outputs/ml_BC_transfer_edge_recovery_v1
```

All hashes recorded in the frozen protocol, proposal, and evaluation manifests
were recomputed after completion and passed. Principal final-output hashes are:

```text
edge_recovery_summary.csv
374770de5aafc1084c874e3fb0b51465b82c1bc015657fc5e2d55d65a9933d95

bisection_cases.csv
bb2669352b896229ade115f6e150534aaffeedfa85f8525d56533e8b43c4d645

ml_BC_edge_recovery.mat
824be85a5ddbc286b2461204ffe63a6a43949d62fbaf5c39038b33ecd4f1219e
```
