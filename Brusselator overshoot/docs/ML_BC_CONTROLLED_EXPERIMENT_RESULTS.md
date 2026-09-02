# Controlled budget-matched B--C boundary experiment

## Research question

Can uncertainty-guided sampling reconstruct the intersection of the B--C basin separator with the frozen two-dimensional landing-aware POD plane more efficiently than a budget-matched space-filling control, and do its proposed brackets recover the independently verified periodic edge orbit (E_{BC})?

This is a controlled first ML experiment, not a full-state basin-learning result. The plane, 25-state seed set, eight-run acquisition budget, 16-state test set, GP family, uncertainty functional, and four edge-bracket proposals per model were frozen before the test outcomes were generated. The exact edge orbit and the held-out C--B--A challenge region were excluded from fitting and acquisition.

## Blind protocol

- Shared training seed: 25 resolved B/C trajectories.
- Added budget: eight active and eight deterministic space-filling trajectories, with one overlapping state (15 unique PDE runs).
- Final model size: 33 labels for each model.
- Fixed test: 16 states, with predictions and model hashes written before any test trajectory was run.
- Dynamical test: four adjacent, label-blind (p_C=0.5) bracket proposals per model, also hash-frozen before endpoint labeling.
- All classifications used the existing physical five-way A/B/C/R5/R6 orbit library, thresholds, and adaptive horizon through at most (t=400).

All 15 acquisition cases, 16 blind tests, and 12 unique bracket endpoints resolved to B or C; none was mapped to `other` or `U`.

## Results

| Quantity | Active | Space-filling control |
|---|---:|---:|
| Blind accuracy | 16/16 | 16/16 |
| Brier score | 0.01594 | 0.02936 |
| Log loss | 0.08682 | 0.12836 |
| Full-plane uncertainty reduction per added trajectory | 0.02385 | 0.01338 |
| Frozen separator-band reduction per added trajectory | 0.002525 | 0.0002803 |
| Verified B/C brackets among four proposals | 0 | 1 |
| First verified bracket rank | -- | 4 |
| Recovered exact (E_{BC}) neighborhood | not tested: no bracket | yes |

The active model reduced the predeclared separator-band uncertainty functional about 9.01 times as much per added trajectory as the control. Its Brier score was about 54% of the control value and its log loss about 68% of the control value. With only 16 test states and identical 16/16 hard-label accuracy, these are descriptive results rather than strong statistical evidence.

The dynamical result reverses the simple ranking. All four active proposals had same-basin endpoint labels. The control model's rank-4 proposal was a genuine B/C bracket. Ten multiclass-safe bisections narrowed it to (9.77\times10^{-4}) of its original adjacent-state segment. The midpoint trajectory then:

- reached a minimum physical phase/reflection-aware distance (1.737\times10^{-4}) from the independently Newton-converged (E_{BC}) orbit;
- remained within the predeclared distance-0.02 edge neighborhood for 175 time units;
- had direct full-state period estimate 3.99772, versus exact period 4.0029149 (relative difference (1.30\times10^{-3}));
- remained unresolved between B and C at the end of the 220-unit diagnostic, as expected for a long edge-shadowing trajectory.

This is strong evidence that the space-filling model produced a dynamically useful bracket recovering the already verified B--C edge mechanism. It is not evidence that its probabilistic separator is globally more accurate, nor that space-filling will generally outperform active learning.

## Interpretation and limitation

The controlled experiment exposes a useful distinction:

> Better classification probabilities and lower integrated posterior uncertainty do not guarantee better proposals for invariant-structure recovery.

The current active acquisition criterion targets entropy, novelty from training states, and batch diversity. It does not model bracket validity, local separator bias, or the probability that a proposed opposite-side pair will produce two different asymptotic PDE outcomes. The dynamical failure is therefore a result to learn from, not a reason to discard the experiment.

Two optimization diagnostics also need follow-up: the active GP signal standard deviation reached its upper bound 30, and the control GP's first length scale reached its upper bound 3. These bounds were predeclared and were not changed after opening the test set. Probability calibration and sensitivity to weak hyperparameter priors or fixed shared kernel parameters should be tested in a new, explicitly versioned experiment.

## Recommended next implementation

The next active-learning iteration should optimize a dynamical-utility acquisition objective rather than entropy alone. Candidate pairs should be ranked using a combination of:

1. opposite-side posterior probability;
2. robustness of that sign change across posterior samples or model ensembles;
3. small state-space separation;
4. novelty relative to existing brackets;
5. expected endpoint-label disagreement;
6. eventual reward when edge tracking approaches a known or newly discovered invariant object.

The present four-proposal result should remain the untouched benchmark. A second preregistered batch can test whether bracket-aware acquisition improves verified-bracket yield while retaining blind proper-score performance.

## Reproducible outputs

- Frozen models and hashes: `experiment_outputs/ml_BC_budget_comparison_v1`
- Blind test evaluation: `experiment_outputs/ml_BC_budget_evaluation_v1`
- Frozen bracket endpoint labels: `experiment_outputs/ml_BC_edge_bracket_labeling_v1`
- Edge recovery: `experiment_outputs/ml_BC_edge_recovery_v1`
- Consolidated table, JSON, and figure: `experiment_outputs/ml_BC_controlled_experiment_v1`

The main machine-readable summaries are `controlled_model_comparison.csv`, `controlled_experiment_summary.json`, `blind_test_individual_predictions.csv`, `frozen_bracket_outcomes.csv`, and `edge_recovery_summary.csv`.
