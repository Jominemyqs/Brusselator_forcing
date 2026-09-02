# Stage 2: phase-aware landing-state geometry

## Scientific purpose

Stage 2 asks how the finite-dimensional set of states deposited when the
forcing returns to `b=10` intersects the state-space organization of the five
currently validated outcomes: stationary A and periodic B, C, R5, and R6.

The projection is exploratory geometry. It is not a learned basin boundary,
an invariant-manifold computation, or evidence that Euclidean proximity alone
determines the final outcome.

The immediate questions are:

1. Does the symmetric-ramp family trace a continuous landing curve?
2. Where do the 16/17 R5--B and 40/45 B--A brackets lie relative to the orbit
   library?
3. Is C geometrically remote from the currently sampled forcing landing set?
4. Do the apparent R6 spatial variants collapse toward one orbit tube or occupy
   visibly different regions?
5. Are there landing states whose reconstruction error or orbit distances make
   them credible novelty candidates?

## Evidence gate inherited from Stage 1

R5 and R6 enter the reference library only because the following finite-time,
`Lx=40` gates have passed:

- direct full-state recurrence after phase-fixed frozen restarts;
- persistence at `N=200,400,600` under tight `ode15s`;
- an `ode45` cross-check at `N=200`;
- local return in two smooth symmetry sectors at amplitudes `1e-4` and `1e-3`;
- agreement of more than one forcing preparation with the same orbit set;
- separation from A, B, C, and from the other recurrence-period class.

This is numerical support for locally attracting periodic outcomes, not a
Floquet proof or a continuum theorem.

## Dataset contract

All states must use `N=400`, `Lx=40`, and the complete concatenated `(u,v)`
field. Every row in the saved manifest must contain:

- source study and raw file;
- protocol identifier and `(bmax,Tup,Thold,Tdown)`;
- forcing-end time and exact landing state;
- final label and label status (`validated`, `provisional`, `novel`, or
  `unresolved`);
- distances to A, B, C, R5, and R6;
- reflection applied during canonicalization;
- phase information when the row belongs to a periodic reference;
- whether the row is one of the 16/17 or 40/45 transition-bracket endpoints.

Pilot and ramp-refinement protocols must be deduplicated by their complete
forcing tuple, not only by their file name. Landing states are taken at the
instant `b(t)` returns to 10, before autonomous relaxation changes the state.

The provisional label must remain available. In particular, ramp 13--15 cases
share the R6 period and move toward the R6 representative but remain slightly
above the strict same-orbit cutoff after the present tight extension. They must
not be silently converted into validated R6 labels.

## Symmetry and phase handling

The Neumann problem has an exact reflection `x -> L-x`; it does not have an
exact continuous translation symmetry. Temporal phase along periodic orbits is
the main nuisance direction.

For a phase-fixed outcome view:

1. Define phase zero by an upward crossing of the orbit-centered spatial mean
   of `v`.
2. If a direct period contains multiple reflection-related crossings, identify
   them modulo reflection and retain one canonical branch.
3. Orient the state with the first nonzero reflection-odd cosine feature. If all
   tested odd features vanish, the state is reflection invariant and no
   orientation is needed.
4. Retain several cycle intersections to display numerical scatter rather than
   presenting a single perfectly sharp point.

For a full-orbit view, orient the complete orbit once using its phase-zero
anchor and apply that same reflection to every phase. Independently reflecting
each phase would introduce artificial discontinuities and could tear a closed
orbit into multiple pieces.

Landing states have no temporal-phase quotient. They receive only the exact
reflection canonicalization.

## POD construction

The primary POD basis is fitted to a balanced validated outcome library using
explicit weights rather than duplicated states. Assign aggregate weight `1/5`
to each of A, B, C, R5, and R6. A receives one sample of weight `1/5`; every
phase sample of a periodic class receives `(1/5)/n_phase`. Subtract the weighted
mean and multiply each centered sample row by the square root of its weight
before SVD. Save all sample and aggregate weights.

Landing states are projected passively and do not determine this first basis.
This prevents a densely sampled forcing family from rotating the coordinates
away from the invariant outcomes being used for interpretation.

Before SVD, concatenate the full fields after applying `sqrt(dx)` quadrature
weighting and separate stored RMS scales for `u` and `v`. Center with the
weighted training mean. Save the mean, scales, right singular vectors,
singular values, explained variance, and all preprocessing metadata.

Report reconstruction errors with 2, 3, 5, 10, and 20 components for every
outcome class and every landing family. If landing-state errors are much larger
than outcome-library errors, fit a secondary landing-aware POD as a diagnostic;
do not replace the primary coordinates silently.

## Required complementary figures

### Phase-fixed representation

Plot the phase-fixed A/B/C/R5/R6 samples, all forcing landing states, final
labels, provisional/novel cases, and emphasized 16/17 and 40/45 brackets. This
is the main figure for comparing outcome neighborhoods and landing geometry.

### Full-orbit representation

Plot uniformly sampled complete periodic orbits as loops, A as a point, and
the same landing states. Produce PC1--PC2 and PC1--PC2--PC3 versions. This view
tests whether the landing set approaches a particular orbit phase, wraps around
an orbit tube, or passes between dynamically distinct neighborhoods.

Also produce an orbit-distance panel and explained-variance/reconstruction-error
panel. A visually separated PCA cluster is not by itself a basin or attractor.

## Quality checks

The Stage 2 builder must fail rather than proceed when:

- a source state is not full `(u,v)` data on `N=400`, `Lx=40`;
- the exact forcing-end state is unavailable;
- duplicate protocols disagree numerically;
- reflecting a state changes its canonical coordinates beyond tolerance;
- an outcome reference fails its saved validation gate;
- a requested label cannot be traced to a documented classifier or validation
  result.

The builder should save a machine-readable state manifest before computing PCA.
This manifest becomes the input contract for later edge-bracket selection and
active learning.

## Decision after Stage 2

Proceed to A--B edge tracking with the 40/45 bracket regardless of whether two
principal components are visually clean, because that bracket is dynamically
validated. Use the 16/17 R5--B bracket second. Introduce nonlinear representation
learning only if POD reconstruction errors or visible folding show that the
linear coordinates are inadequate; class overlap alone is not sufficient.

## Reflection-continuity diagnostic

Feature-based reflection canonicalization can jump when a reflection-odd
feature crosses zero. Therefore every landing row must retain four versions:

1. the raw forcing-generated state;
2. its exact reflected copy;
3. the protocol-independent feature-canonical state;
4. a continuity-oriented state for ordered one-parameter families.

For a family ordered by its varied parameter, choose the raw or reflected state
that minimizes physical full-state distance to the previously oriented family
member. Store every orientation decision and plot the raw, feature-canonical,
and continuity-oriented curves separately. The first Stage 2 test is whether
the raw forcing-to-landing map is continuous and how reflection quotienting
represents it. A quotient-coordinate jump must not be interpreted as a fold or
disconnected branch.

## Phase-section validation

For each periodic reference, the builder must record every upward crossing of
the orbit-centered spatial mean of `v` over one direct period and compute

`d/dt (mean(v)-orbit_mean(mean(v)))`

at each crossing. The report must test transversality, crossing-count stability,
temporal-resampling stability, and branch reproducibility. Multiple crossings
must be compared in the physical phase/reflection-aware norm. If they are not
reflection-equivalent, retain all branches or fall back to closest-template
phase minimization. A branch may not be selected silently.

## Separation of norms

Three constructions remain distinct and are saved under different names:

- the quadrature-weighted physical full-state norm in original `u,v` units;
- the physical norm minimized over phase/reflection for classification;
- the component-scaled norm used only for POD and reconstruction diagnostics.

No scientific orbit-distance claim may be made from POD coordinates. Alongside
the primary separately scaled POD, compute an alternative POD using the same
quadrature weights but no separate `u/v` normalization. Report subspace or
coordinate sensitivity so visual conclusions are not artifacts of scaling.

## Implemented result

The specification is implemented by `build_stage2_landing_geometry.m` and the
finalized output is `experiment_outputs/stage2_landing_geometry_v2`. The
versioned result includes the state and family manifests, phase-section audit,
balanced POD weights, unscaled-component sensitivity, physical orbit-distance
matrix, reconstruction diagnostics, and all required figures. The declared
landing-error fallback was triggered, so the output also contains a secondary
50/50 outcome/landing-weighted POD diagnostic. It is labeled diagnostic-only
and does not replace the balanced outcome-only primary coordinates.
