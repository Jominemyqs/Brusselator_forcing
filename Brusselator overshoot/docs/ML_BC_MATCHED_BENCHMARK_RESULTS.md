# Repeated matched benchmark for task-aware B--C edge-bracket acquisition

> **Audit status (2026-08-15).** The bracket-yield comparison in this preserved
> matched-v1 record is superseded for manuscript inference. The frozen proposal
> script used the latent GP mean as if it were `p(C)` when computing pointwise
> entropy and local orientation. The PDE label bank, joint opposite-sign
> probability, guarded refinements, and exact `E_BC` recovery are unaffected.
> See `ML_BC_MATCHED_BENCHMARK_PROBABILITY_INDEX_AUDIT.md` and the corrected
> `ML_BC_ACQUISITION_ROBUSTNESS_RESULTS.md`. The numerical material below is
> retained unchanged as the historical versioned record.

## Status and frozen scope

The versioned benchmark `ml_BC_matched_dynamical_benchmark_v1` completed on
2026-08-14. Its protocol, matched initial designs, acquisition rules, endpoint
budgets, proposal batches, and success criteria were frozen before candidate
outcomes were opened for evaluation.

The experiment addresses one controlled methodological question:

> Given sparse labels on a fixed two-dimensional B--C initial-state family, does
> a pair-aware probabilistic acquisition rule produce usable edge-tracking
> brackets more reliably than pointwise uncertainty or geometric space filling
> under the same endpoint-label budget?

It is a robustness benchmark for bracket acquisition on the known
`N=400`, `Lx=40`, frozen-`b=10` B--C development plane. It is not a full-state
basin reconstruction, an independent physical experiment, a prospective
discovery of the B--C mechanism, or a test of general rare-basin discovery.

## Dynamical motivation

Ordinary active classification assigns utility to individual states. Edge
tracking instead requires a relation between two states: their asymptotic
outcomes must differ, they must be sufficiently close to define a useful local
bracket, and guarded refinement must not expose an intervening third outcome.
The dynamical task therefore induces the learning objective:

```text
sparse PDE labels
    -> nearby opposite-outcome pair
    -> multiclass-safe bracket
    -> edge tracking
    -> independently verified invariant object.
```

The benchmark compares acquisition methods at the first two arrows and then
uses the known periodic B--C edge orbit `E_BC` as an independent dynamical
target for conditional verification.

## Frozen experimental design

- State family: the existing `17 x 17` landing-aware PC1--PC2 plane.
- Frozen system: `b=10`, `N=400`, `Lx=40`.
- Matched initial-design replicates: 12.
- Training labels per replicate: 25.
- Fixed anchors: state 1, labeled B, and state 289, labeled C.
- Other training states: 23 outcome-blind randomized-maximin selections.
- Acquisition methods:
  - `space_filling`;
  - `pointwise_uncertainty`;
  - `pair_aware`.
- Acquisition batch: four disjoint pairs per method and replicate.
- New endpoint-label budget: eight labels per method and replicate.
- Primary independent comparison unit: the matched initial-design replicate,
  not an individual pair proposal.

Within each replicate, all methods receive the same initial training design,
candidate coordinates, admissible short-pair set, GP fit, pair count, endpoint
budget, and midpoint-diversity rule.

The classifier is the frozen ARD squared-exponential probit GP with Laplace
inference. `pointwise_uncertainty` uses the geometric mean of the two endpoint
entropies together with the same short-distance and local-orientation factors
used by the pair-aware method. `pair_aware` replaces marginal entropy with the
joint posterior probability that the two latent GP values have opposite signs.

## Label bank and blinding audit

The complete fixed-plane bank contains 289 resolved PDE outcomes:

| Final outcome | Count |
|---|---:|
| B | 159 |
| C | 130 |
| other | 0 |
| unresolved | 0 |

The binary B/C plane gate therefore passed for this controlled family. Of the
289 labels, 76 were preserved from prior versioned experiments and 213 missing
states were generated for this benchmark. The recorded sum of unique PDE
labeling runtimes is `13529.18 s`.

The full bank was used as a masked replay resource. Although its labels had
been precomputed, the acquisition implementation did not expose candidate
outcomes. The frozen proposal manifest records:

- proposal batches frozen: `true`;
- candidate outcomes used: `false`;
- only the current replicate's training labels exposed: `true`;
- exact edge orbit used: `false`;
- held-out C--B--A region used: `false`.

The proposals were hashed before endpoint outcomes were opened. Consequently,
the repeated experiment measures label efficiency under hidden-label replay.
It estimates how an online strategy would allocate its label budget on this
fixed bank, but it is not a realized end-to-end wall-clock saving: the PDE bank
was computed in full to make repeated matched comparisons possible.

## Primary result: bracket discovery

An acquisition proposal is a valid bracket when its two independently opened
endpoint labels are B and C in either order.

| Method | Replicates | Proposals | Endpoint-label budget | Valid B--C brackets | Pooled bracket yield | Mean brackets per replicate | Replicates with a bracket | Brackets per endpoint label |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| Pair-aware | 12 | 48 | 96 | 32 | 0.6667 | 2.6667 | 12/12 | 0.3333 |
| Pointwise uncertainty | 12 | 48 | 96 | 16 | 0.3333 | 1.3333 | 10/12 | 0.1667 |
| Space filling | 12 | 48 | 96 | 4 | 0.0833 | 0.3333 | 4/12 | 0.0417 |

Thus, on this benchmark, pair-aware acquisition doubles pooled bracket yield
relative to pointwise uncertainty and increases it eightfold relative to space
filling. Because proposals within a four-pair batch are correlated, these
pooled ratios are descriptive. The paired initial-design comparison is the
primary analysis.

### Paired replicate-level comparisons

| Comparison | Mean paired difference in valid brackets | Paired bootstrap 95% interval | Win/tie/loss for pair-aware | Exact two-sided sign-flip `p` |
|---|---:|---:|---:|---:|
| Pair-aware minus pointwise uncertainty | 1.3333 | [0.8333, 1.8333] | 9/3/0 | 0.00390625 |
| Pair-aware minus space filling | 2.3333 | [1.6667, 3.0000] | 11/1/0 | 0.0009765625 |

Pair-aware acquisition never loses to either comparator across the 12 matched
initial designs. The repeated result therefore replaces the original
four-proposal ordering as the primary evidence for a task-aware acquisition
advantage.

### Labels to the first bracket

Pair-aware acquisition finds a first bracket after a median of two revealed
endpoint labels. Among the replicates in which the comparator succeeds, the
corresponding conditional medians are three labels for pointwise uncertainty
and seven for space filling. The comparator medians must remain explicitly
conditional because pointwise fails to find any bracket in 2/12 replicates and
space filling fails in 8/12.

## Conditional guarded refinement and edge recovery

The downstream test is deliberately conditional on bracket discovery. For each
method and replicate that produces at least one valid bracket, the frozen rule
selects its highest-ranked valid proposal. This gives:

| Method | Conditional selections | Multiclass-safe refinements | Exact `E_BC` recoveries |
|---|---:|---:|---:|
| Pair-aware | 12 | 12 | 12 |
| Pointwise uncertainty | 10 | 10 | 10 |
| Space filling | 4 | 4 | 4 |
| Total method--replicate selections | 26 | 26 | 26 |

Six selections duplicate an identical physical endpoint pair across methods or
replicates. Exact deduplication leaves 20 unique physical pairs, each integrated
only once and attributed back to every frozen selecting group.

Every unique pair:

- completes all 10 prescribed bisection steps;
- retains a B endpoint and a C endpoint;
- produces only B or C midpoint classifications;
- reaches final bracket width `0.0009765625` in interpolation coordinate;
- enters the independently defined neighborhood of the Newton-converged
  periodic edge orbit `E_BC`.

Across the 20 unique recoveries:

- the minimum phase/reflection-aware physical orbit distance ranges from
  `4.9903e-5` to `2.9317e-4`, with median `1.0953e-4`;
- the edge-shadow duration ranges from `138.5` to `186.5`, with median `170.5`;
- the direct period estimate has relative error between `1.4852e-4` and
  `1.9236e-3` relative to
  `T_E = 4.0029149135788`.

The candidate final outcome is recorded as `unresolved` in all 20 cases. This
is expected: the candidate trajectories shadow an unstable edge orbit for a
long interval rather than rapidly converging to B or C. The recovery claim is
based on the independently prescribed orbit-distance, shadow-duration, and
direct-period diagnostics, not on assigning the edge trajectory an attractor
label.

The recorded runtime of the 20 deduplicated guarded-refinement and candidate
recoveries is `18298.57 s`. Together with the recorded label-bank runtime, the
sum is `31827.74 s`; this is accumulated simulation runtime, not necessarily
elapsed wall time.

## Scientific interpretation

The repeated benchmark supports three statements.

1. **The acquisition objective matters.** On the fixed B/C plane and under the
   same revealed-label budget, optimizing the joint event needed by edge
   tracking produces valid brackets more efficiently and consistently than
   optimizing endpoint uncertainty or geometric novelty.
2. **Pointwise uncertainty is useful but not task-optimal here.** It produces
   valid brackets in 10/12 replicates and should not be characterized as a
   general failure. Its bracket yield is nevertheless half that of pair-aware
   acquisition.
3. **Dynamical verification remains separate from ML.** Conditional on finding
   a valid bracket, all three methods recover the same independently computed
   periodic orbit. The ML result concerns access to a useful initialization;
   the PDE integration, guarded refinement, Newton solve, and Floquet analysis
   establish the invariant mechanism.

The method-specific conditional orbit-distance and shadow-duration medians are
not an unconditional ranking of acquisition methods. Each method contributes a
different selected subset, and the primary methodological comparison is valid
bracket yield at matched budgets.

## Relationship to the original pilot and held-out challenge

The original four-proposal pilot produced:

| Method | Valid brackets |
|---|---:|
| Pointwise uncertainty | 0/4 |
| Space filling | 1/4 |
| Pair-aware | 3/4 |

That result remains useful as the observation that motivated the pair-aware
objective, but it is too small to support a comparative claim by itself. The
repeated experiment shows that the pilot's complete pointwise failure was not
representative: pointwise uncertainty is intermediate between pair-aware
acquisition and space filling across varied initial designs.

The frozen eight-run C--B--A held-out experiment remains unchanged:

```text
7 A + 1 C + 0 B.
```

The learner localized the relevant region but did not directly sample the
narrow B intrusion under the prespecified budget. Post-budget guarded
refinement subsequently found B and recovered the same edge orbit, but that
diagnosis is not counted as prospective success. This negative result is the
necessary limit on the controlled benchmark: efficient B/C bracket discovery
on a known binary plane does not establish rare-basin or novelty discovery in
multibasin geometry.

## Frozen manuscript claims

The manuscript may state:

> In a controlled repeated benchmark on the fixed B--C development plane,
> task-aware pair acquisition doubles valid-bracket yield relative to pointwise
> uncertainty and increases it eightfold relative to space filling under equal
> endpoint-label budgets. It produces at least one bracket in every matched
> initial design and never loses to either comparator at the replicate level.

It may also state:

> Every conditionally selected valid bracket passes multiclass-safe refinement
> and recovers the independently Newton-converged periodic edge orbit, showing
> that the acquisition proposals lead to the same verified invariant mechanism.

The manuscript must not state:

- that the 12 replicates are independent physical experiments;
- that the benchmark demonstrates realized wall-clock savings;
- that pointwise uncertainty cannot find edge brackets;
- that pair-aware acquisition changes the edge state reached once a valid
  bracket is supplied;
- that the learner discovers `E_BC` without classical edge tracking;
- that the method reconstructs the full basin boundary;
- that it reliably discovers narrow rare basins or C-accessible forcing
  protocols;
- that the result generalizes beyond the fixed two-dimensional `N=400`,
  `Lx=40` development plane without further experiments.

## Manuscript placement and frozen figure content

This result belongs after the exact B--C edge state has been established. The
section should open with the edge-tracking requirement, which creates the
pair-acquisition problem, rather than with a generic introduction to ML.

The main benchmark figure should contain:

1. a schematic contrasting isolated-point uncertainty with pairwise bracket
   utility;
2. replicate-level matched bracket counts, with the same initial design linked
   across methods;
3. bracket yield per proposal and per endpoint label;
4. the conditional funnel from valid proposal to guarded B/C refinement to
   exact `E_BC` recovery.

The original pilot belongs in an inset or supplement. The held-out C--B--A
failure should remain a separate panel or figure so that controlled success is
not conflated with prospective rare-basin discovery.

## Reproducibility and integrity records

Critical verified files and SHA-256 hashes are:

| Artifact | SHA-256 |
|---|---|
| Complete 289-state label bank | `3f1e9fef80b33923faf45b844446b693ee11f5d9fac8af9554cd8d547eebe0c5` |
| Label-source audit | `b9b288721ded4c0e99e718b2eb9cc67f7ae86b380280349243d9b314ff9d37c2` |
| Frozen protocol manifest | `12497001f031fbe73c0f3b939be5bd1872dd972bdf9939a5fffd0e1b75bf59b7` |
| Frozen shard plan | `499168d80a201ef777caa830eb960b261e264acc7b09ba61694fdf300732a8c6` |
| Matched initial designs | `68703376002c230b9b3c24537e851c6c109bb02b72f3cd00ca7a9298d1c8c964` |
| Frozen proposal manifest | `ed1b22be17f186dc97e05ddd6766ef53f48efe9723b95b3bb6ed629d33cdf111` |
| Frozen pair proposals | `202ac3bc5ce841b7e9fd0febf323fd4172469fcec8273e0e55e6fbcf25e0a94e` |
| Selected conditional pairs | `cc97f7204e36a2865cb9a7650b374b0fd8b574741f85f51fb820c45bc852ee98` |
| Unique physical pairs | `2e0c163bcf534f2c7c2b93f87465af7537fe7b7c52230867c2d492c0465a618c` |
| Edge-recovery attribution | `f1e3a39095366b03a11759855a592317d157e14aa0116bb7d07b851fc0ed035d` |

The completion audit independently verified all listed hashes. It also verified:

- 289 complete bank rows with exactly 159 B and 130 C outcomes;
- 144 frozen proposals: 48 per method over 12 replicates;
- 200 bisection rows: 10 for each of 20 unique physical pairs;
- exactly 100 B and 100 C bisection midpoint outcomes;
- no remaining partial files and no missing candidate or bisection outputs;
- 26/26 attributed selections marked multiclass-safe and exact-edge recovered.

The repository is not currently a Git repository, so the generated metadata
correctly records `git_available=false` and no commit identifier.

## Machine-readable outputs

- Protocol and matched designs:
  `experiment_outputs/ml_BC_matched_benchmark_protocol_v1`
- Complete PDE label bank:
  `experiment_outputs/ml_BC_matched_label_bank_v1`
- Frozen outcome-blind proposals:
  `experiment_outputs/ml_BC_matched_acquisition_proposals_v1`
- Opened endpoint outcomes and replicate-level comparison:
  `experiment_outputs/ml_BC_matched_acquisition_evaluation_v1`
- Deduplicated guarded refinement and edge diagnostics:
  `experiment_outputs/ml_BC_matched_edge_recovery_v1`
- Final conditional dynamical summary:
  `experiment_outputs/ml_BC_matched_dynamical_benchmark_v1`

The permanent top-level machine-readable result is
`experiment_outputs/ml_BC_matched_dynamical_benchmark_v1/matched_dynamical_benchmark_summary.json`.
