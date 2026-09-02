# Initial probabilistic B--C boundary model

## Question and scope

`ml_BC_probabilistic_boundary_v1` asks whether a probabilistic learner trained
only on the 25 frozen-PDE seed labels places high uncertainty along the observed
B--C separator intersection and can propose new, physically simulated states
near that region.

This is a two-dimensional initial-state experiment. It is not a learned PDE
solver, a full-state basin classifier, or a computation of the stable manifold.
The exact periodic edge orbit and held-out C--B--A challenge region are excluded
from fitting and acquisition.

## Model

The classifier is an auditable local implementation of a binary Gaussian
process with:

- normalized `(alpha_1, alpha_2)` coordinates;
- an anisotropic squared-exponential kernel;
- a probit observation model;
- Laplace posterior inference;
- marginal-likelihood hyperparameter fitting from four deterministic starts;
- no random seed or random acquisition tie-breaking.

The optimized full-data parameters are

```text
lengthscale_1 = 2.67106
lengthscale_2 = 1.36708
signal standard deviation = 16.2435
log marginal likelihood = -7.35968
```

All parameters are interior to their declared optimization bounds. A preliminary
cap of 10 was rejected before freezing because the signal scale reached that
bound; increasing the declared upper bound to 30 gives a stable interior optimum
near 16.24.

The fitted posterior separates all 25 training labels. Nested leave-one-out
refitting gives

```text
accuracy = 0.96
Brier score = 0.05824
log loss = 0.22342
```

The only leave-one-out error is state 285, the B point at
`(alpha_1,alpha_2)=(0.0406483,0.0352625)` that exposed the strongest coarse
separator curvature. Its held-out probability is `p_C = 0.56465`.

These are internal small-sample diagnostics. The saved probabilities are
Laplace-probit posterior probabilities, but empirical calibration has not been
established. The 16 fixed-test states remain unlabeled; their blind predictions
have been frozen before their outcomes are computed.

## Boundary and acquisition

The dense posterior has one connected `p_C = 0.5` contour. Its geometry bends
toward positive `alpha_2` as `alpha_1` increases, matching the qualitative
curvature visible in the 5-by-5 PDE labels.

The maximum acquisition-pool entropy is `0.999966` bits. A reproducible
eight-state batched active proposal combines posterior entropy with distance
from the training set and diversity within the proposed batch:

```text
111, 24, 251, 198, 287, 164, 58, 42
```

The first three have probabilities `p_C = 0.5035`, `0.5072`, and `0.4835`.
The remaining points add spatial coverage along the uncertain band rather than
duplicating nearly identical queries.

The matched eight-state comparison batch uses deterministic farthest-point
space filling without classifier uncertainty:

```text
20, 24, 28, 32, 88, 92, 96, 100
```

State 24 is selected by both methods, so the two eight-trajectory budgets share
one physical simulation and contain 15 unique new states. This overlap is
retained rather than altered after seeing the proposals.

## Reproducibility and next test

The model implementation, synthetic numerical test, run configuration, fitted
arrays, optimization records, nested leave-one-out predictions, full candidate
probabilities, blind test predictions, contour, and acquisition lists are all
saved. Twenty post-fit integrity checks pass.

The next controlled computation is to label the active and space-filling batches
with the same adaptive frozen-PDE runner. Two models can then be refit from the
same 25-state seed plus their respective eight-state additions. The fixed test
labels should remain unexamined until those comparison models and evaluation
metrics are frozen.
