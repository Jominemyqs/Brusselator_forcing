# ML-assisted B--C boundary discovery: first experiment

## Research question

Can uncertainty-guided learning locate the intersection of a known basin
separator with a physically motivated low-dimensional state family, and do
ML-selected boundary brackets dynamically converge to the independently
verified periodic edge state?

The purpose is not to learn the Brusselator time stepper. The PDE remains the
ground-truth simulator. ML proposes informative initial states and possible
boundary brackets; edge tracking verifies the invariant dynamics.

## Controlled validation problem

Use the clean B--C region first. The exact B and C attractors, periodic edge
orbit, Floquet spectrum, unstable direction, and branch outcomes have now been
verified independently and across `N = 200, 400, 600`. None of this invariant
information is supplied as a training label beyond the eventual B/C outcome of
each sampled initial state.

Construct a two-dimensional family from the existing landing-aware Stage-2 POD
basis,

```text
X(alpha_1,alpha_2) = X_0 + alpha_1 phi_1 + alpha_2 phi_2.
```

`X_0` should be a forcing-landing or bracket-midpoint state in the clean B--C
region, not the exact edge orbit. The exact edge orbit is reserved for
evaluation. POD scaling is used only to define the sampling plane; all outcome
and edge distances remain in the documented physical quadrature-weighted norm.

Save every full initial state and reject or flag samples that violate declared
physical admissibility checks. Reflection canonicalization must retain the raw,
reflected, and chosen orientations so quotient discontinuities remain auditable.

## Labels

For every selected point, integrate the actual frozen PDE and assign

```text
B       approaches the exact B orbit,
C       approaches the exact C orbit,
other   approaches another validated outcome,
U       unresolved at the adaptive horizon.
```

The classifier must never convert `other` or `U` into the nearest B/C label.
Initially train the probabilistic B/C boundary model only on resolved B/C
samples while retaining a separate novelty/unresolved indicator.

## First model

Use the simplest calibrated probabilistic model that works at the observed
sample size, preferably a Gaussian-process classifier in the two POD
coordinates. Do not begin with an autoencoder, Deep Koopman model, or learned
PDE surrogate.

The primary object is the predictive uncertainty field, not overall accuracy:

```text
p_B(z), p_C(z), predictive entropy, and the p_B approximately p_C contour.
```

Compare a fixed initial design with sequential entropy-guided sampling. Use the
same initial seed set and simulation budget when comparing active and uniform
sampling.

## Dynamical evaluation

The main evaluation is not classification percentage. For ML-selected uncertain
regions:

1. find nearby resolved B/C points with opposite labels;
2. use their complete PDE states as edge-tracking brackets;
3. dynamically rebracket with the established physical classifier;
4. measure phase/reflection-aware distance from the recovered edge trajectory
   to the exact period-`4.003` B--C edge orbit.

Report:

- PDE trajectories required before the separator is detected;
- distance between the learned uncertainty contour and a dynamically refined
  reference boundary in the same plane;
- fraction of ML-proposed brackets that remain valid B/C brackets;
- fraction that recover the verified periodic edge orbit;
- edge-orbit distance and shadowing duration for each recovered trajectory;
- active-versus-uniform efficiency at matched simulation budgets;
- frequency and location of `other` and `U` outcomes.

The scientifically decisive validation question is:

> Does classifier uncertainty trace the stable-manifold intersection of the
> verified periodic B--C edge state?

## Held-out discovery problem

Do not use the narrow interval near `lambda = 0.276` to design or tune the
initial method. After the clean B--C validation is frozen, apply the same
uncertainty/novelty strategy to the nominal C--A region.

The held-out question is whether the method independently discovers the narrow
B-labeled region already indicated by the coarse guarded refinement:

```text
C -> B -> A.
```

Success here means detecting an unexpected third-basin intrusion and proposing
new simulations that resolve it. It does not require or justify calling the
structure fractal, riddled, or a direct C--A boundary.

## Implementation gates

1. Freeze the two-dimensional sampling plane, center, scaling, bounds, and
   physical admissibility rules before outcome simulations.
2. Produce a small uniform reference design and explicit train/test splits.
3. Calibrate the probabilistic B/C model and retain `other/U` separately.
4. Implement budget-matched uniform and uncertainty-guided acquisition.
5. Validate selected brackets through dynamic edge tracking and exact-edge
   distance, not only prediction accuracy.
6. Freeze the method before exposing the held-out `lambda approximately 0.276`
   challenge region.

Only if the two-dimensional POD plane visibly folds distinct dynamical regions
onto one another should nonlinear representation learning be considered.

## Implementation status: frozen geometry complete

The first two implementation gates are now instantiated by
`ml_BC_boundary_geometry_v1`:

- `X_0` is the reflection-aligned midpoint of the independently classified
  `lambda = 0.1875` B and `lambda = 0.187890625` C bracket states.
- The physical sampling directions are decoded from landing-aware POD
  components 1 and 2.
- Coordinate bounds are fixed from the observed score ranges of validated
  forcing-generated B landings, followed by an explicit positivity gate.
- A 17-by-17 master grid fixes a 25-point initial seed, 16-point untouched test
  set, and 248-point acquisition pool.
- The exact edge orbit and held-out `lambda approximately 0.276` region are
  excluded from all construction sources.
- No generated state has yet been simulated or labeled.

The known B--C bracket displacement has `0.975878` of its scaled norm in the
selected two-mode plane. Reflection canonicalization also produces a visible
coordinate switch, so raw, reflected, canonical, and continuity-oriented states
are all retained. Numerical details and limitations are recorded in
`ML_BC_BOUNDARY_GEOMETRY_RESULTS.md`.

## Implementation status: initial adaptive labels complete

`ml_BC_boundary_initial_labeling_v1` has now evolved all 25 common initial-seed
states with the frozen PDE. It uses exact Newton-converged B/C references and
validated A/R5/R6 references, with the original physical classifier thresholds.
Unresolved cases are extended from 80 in 40-unit chunks to at most 400.

The seed contains 14 B and 11 C labels, with no `other` or `U` outcomes. No case
resolved at the initial 80-unit horizon; the longest case required 280 units.
The coarse 5-by-5 map contains seven nearest-neighbor B/C brackets and shows a
curved separator intersection. The fixed 16-point test set and all 248
acquisition candidates remain unlabeled. See
`ML_BC_BOUNDARY_LABELING_RESULTS.md`.

## Implementation status: first probabilistic model complete

`ml_BC_probabilistic_boundary_v1` fits an ARD squared-exponential probit GP by
Laplace inference using only the 25 B/C seed labels. Deterministic nested
leave-one-out refitting gives 96% accuracy; the single error is the B point that
revealed the strongest coarse boundary curvature. The posterior has one curved
`p_C = 0.5` contour and reaches almost one bit of entropy on the acquisition
pool.

An eight-state entropy/novelty/diversity batch and an eight-state deterministic
space-filling comparison batch are now frozen. Their budgets share one state,
giving 15 unique proposed PDE runs. Blind probabilities for the 16 fixed-test
points are saved, but their labels remain uncomputed. Exact-edge and held-out
challenge data remain excluded. See
`ML_BC_PROBABILISTIC_BOUNDARY_RESULTS.md`.
