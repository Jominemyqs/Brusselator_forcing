# Frozen robustness design for B--C edge-bracket acquisition

## Scientific scope

This versioned study is a retrospective masked replay on the existing complete
`17 x 17` B--C PDE-label bank. It tests a deliberately conditional question:

> Given a prelocalized state family containing a candidate pairwise separator,
> does a task-aware acquisition objective generate usable opposite-outcome
> edge-tracking brackets more efficiently than pointwise uncertainty or space
> filling, and is that ordering stable under bounded benchmark perturbations?

It does not test global basin discovery, rare-basin discovery, transfer to an
independent PDE state family, or realized end-to-end wall-clock savings. No new
PDE integration is authorized by this protocol.

## Frozen matched designs

- Forty randomized-maximin designs are generated from geometry alone.
- State 1 (B) and state 289 (C) are fixed anchors.
- Each design contains an ordered sequence of 35 states.
- The `n0=15`, `n0=25`, and `n0=35` designs are nested prefixes of that same
  sequence, which makes initial-budget comparisons matched within replicate.
- Candidate outcomes remain masked during proposal construction.

## Reference experiment

The reference condition retains the primary benchmark choices:

- initial labels: `n0=25`;
- kernel: ARD squared exponential;
- pair score: `qDO`;
- distance scale: `0.18`;
- admissible normalized pair separation: `[0.1249999, 0.2500001]`;
- four disjoint proposed pairs and eight newly revealed endpoint labels.

Here `q` is the joint posterior probability that the two latent GP values have
opposite signs, `D` is the short-pair distance penalty, and `O` is the
local-normal orientation factor.

## One-factor-at-a-time perturbations

The reference is perturbed along four prespecified axes:

1. Initial labels: `n0=15,25,35`.
2. Kernel: squared exponential, Matérn `3/2`, and Matérn `5/2`.
3. Pair-score ablation: `q`, `qD`, and `qDO`.
4. Soft distance scale: `0.75`, `1.0`, and `1.25` times the reference value.

The hard admissible pair-distance interval is fixed in normalized POD-plane
coordinates. Thus the distance experiment changes only the declared soft
penalty and does not redefine a grid cell as a physical distance.

This is not a full factorial study. The one-factor-at-a-time construction keeps
the interpretation bounded and avoids turning a single physical plane into an
exercise in post hoc hyperparameter optimization.

## Frozen comparison and outputs

Every scenario compares pair-aware acquisition with pointwise uncertainty and
space filling under the same matched design, GP fit, candidate family, four-pair
batch, and eight-label endpoint budget. The initial-design replicate is the
independent analysis unit.

Primary outputs are:

- valid B--C brackets among four proposals;
- probability of at least one valid bracket by 2, 4, 6, and 8 revealed labels;
- brackets per revealed endpoint label;
- paired wins, ties, and losses;
- paired bootstrap intervals and exact sign tests;
- the `q`, `qD`, `qDO` ablation ordering.

## Probability-index audit and correction

Before this extension was frozen, an audit found that the original matched-v1
proposal script assigned `predict_proba(...)[1]` to its probability field. In
the preserved GP interface, index 0 is the class probability and index 1 is the
latent mean. The saved v1 table confirms the mismatch: 96 of 144 proposal rows
contain at least one reported “probability” outside `[0,1]`.

This affects the v1 pointwise entropy comparator and the probability-gradient
orientation factor used by pointwise and pair-aware scoring. It does not affect
the PDE label bank, the joint latent opposite-sign probability `q`, geometric
space filling, guarded edge tracking, or the exact periodic edge orbit. The
present protocol corrects the index and therefore supersedes the v1
pair-aware-versus-pointwise bracket-yield comparison for manuscript inference.
The original versioned artifacts remain unchanged as an audit record.

This robustness suite does not repeat guarded edge tracking for every replayed
proposal. The primary benchmark already established that all 20 unique,
conditionally selected valid brackets passed guarded refinement and recovered
the independently verified periodic orbit `E_BC`. The extension tests the
upstream acquisition claim; it does not manufacture additional dynamical
replicates from identical saved PDE states.

## Interpretation gate

Even if pair-aware acquisition is favored under every perturbation, the maximum
allowed claim remains:

> Within this prelocalized B--C state family, task-aware pair acquisition is
> robustly more efficient at producing dynamically usable bracket candidates
> than pointwise uncertainty or space filling under matched label budgets.

The prospective held-out C--B--A failure remains unchanged and is the evidence
that this conditional result does not automatically extend to narrow unseen
basin components.
