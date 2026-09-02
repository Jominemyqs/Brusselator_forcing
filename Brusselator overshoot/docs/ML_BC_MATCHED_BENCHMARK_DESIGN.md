# Repeated matched benchmark for edge-bracket acquisition

## Dynamical motivation

Edge tracking does not require an isolated uncertain state. It requires two
nearby states with different asymptotic outcomes, followed by a guarded check
that their connecting segment is not obstructed by another known basin. The
learning problem is therefore induced by the dynamical task:

`multibasin geometry -> opposite-outcome pairs -> edge brackets -> invariant-object recovery`.

The original four-proposal comparison is retained unchanged as a proof of
concept. This versioned experiment tests whether its ordering is robust to the
sparse initial design.

## Frozen repeated design

- Controlled state family: the existing 17-by-17, two-dimensional, landing-aware
  POD plane at `N=400`, `Lx=40`, and frozen `b=10`.
- Replicates: 12 matched initial designs.
- Training budget: 25 labels per replicate.
- Fixed anchors: state 1 (`B`) and state 289 (`C`).
- Remaining design points: 23 outcome-blind randomized-maximin selections,
  using a saved seed and uniform choice among the 12 most distant available
  states at each step.
- Acquisition methods: geometric space filling, marginal pointwise
  uncertainty adapted to pair proposals, and joint-posterior pair acquisition.
- Acquisition budget: four disjoint pairs, hence eight newly revealed endpoint
  labels, per method and replicate.
- Matched comparison: all methods receive the same initial design, candidate
  coordinates, admissible pair set, GP fit, and pair-count budget within each
  replicate.

The complete design, configuration, hashes, and the 213-state missing-label
manifest were frozen before those remaining outcomes were generated. The 76
previously available labels are preserved rather than recomputed.

## Acquisition definitions

All proposed pairs have normalized separation between approximately one and
two grid spacings and all four pairs in a batch have disjoint endpoints.

- `space_filling` ranks pairs using endpoint novelty relative to the training
  design and short-pair geometry.
- `pointwise_uncertainty` ranks pairs using the geometric mean of their two
  marginal predictive entropies, with the same short-pair and local-normal
  factors used by the pair-aware rule.
- `pair_aware` replaces marginal entropy with the joint posterior probability
  that the two latent GP values have opposite signs.

All three use the same midpoint-diversity rule when constructing the four-pair
batch. The GP is an ARD squared-exponential probit classifier with Laplace
inference and a frozen weakly regularized hyperparameter fit.

## Evaluation hierarchy

The initial-design replicate is the primary independent unit. Individual
brackets from one batch are not treated as independent statistical replicates.

1. Endpoint bracket validity: one endpoint reaches B and the other C.
2. Bracket yield per four proposals and per eight revealed PDE labels.
3. Time to the first valid pair in two-label increments.
4. Ten-step multiclass-safe refinement of the highest-ranked valid pair for
   every method/replicate.
5. Conditional recovery of the independently Newton-converged periodic edge
   orbit `E_BC`.

Paired method differences are summarized across the 12 matched designs using
descriptive intervals, win/tie/loss counts, and an exact sign-flip test. These
algorithmic replicates share one precomputed PDE-label bank and therefore do
not constitute 12 independent physical experiments.

## Interpretation limits

This is a controlled robustness benchmark on the known B--C development plane.
It evaluates whether task-aware acquisition is more reliable for generating
edge-tracking initializations. It is not a full-state basin reconstruction, an
independent discovery of the B--C mechanism, or evidence that the method can
generally discover rare third-basin intrusions. The frozen negative C--B--A
challenge remains the prospective limitation test.

Machine-readable protocol files are in
`experiment_outputs/ml_BC_matched_benchmark_protocol_v1`.
