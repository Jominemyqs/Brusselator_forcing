# Frozen two-dimensional ML boundary geometry

## Purpose

`ml_BC_boundary_geometry_v1` completes the pre-simulation gate for the first
ML-assisted basin-boundary experiment. It defines a reproducible family of
initial states at frozen `b = 10`, `Lx = 40`, and `N = 400` without using the
exact periodic edge orbit and without assigning outcome labels to newly
generated states.

The object constructed here is a sampling plane, not a learned basin boundary:

```text
X(alpha_1, alpha_2) = X_0 + alpha_1 phi_1 + alpha_2 phi_2.
```

## Center and directions

The center `X_0` is the complete-state midpoint of the independently classified
B/C bracket at

```text
lambda_B = 0.1875,
lambda_C = 0.187890625.
```

The C endpoint is first compared with its exact spatial reflection and aligned
to the B endpoint; the midpoint is then feature-canonicalized once. The exact
edge orbit is not loaded by the builder.

`phi_1` and `phi_2` are the first two modes of the saved Stage-2 landing-aware
POD. They are decoded from the stored scaled, trapezoidal coordinates into
physical `u` and `v` perturbations. The two modes explain `57.2111%` and
`12.1882%` of that POD training covariance, respectively (`69.3993%` total).
This percentage describes the POD training measure; it is not a claim that the
full basin geometry is two-dimensional.

The known bracket displacement is unusually well aligned with this plane. Its
two-mode scaled projection retains `0.975878` of the displacement norm. The
projected endpoint coordinates relative to the midpoint are

| Endpoint | alpha_1 | alpha_2 |
|---|---:|---:|
| B | `-1.96454e-5` | `-2.53396e-4` |
| C | `+1.96454e-5` | `+2.53396e-4` |

The saved B/C bracket separation is `2.11478e-4` in the symmetric physical
relative full-state norm. Endpoint labels validate the source bracket only;
they are not inherited by their two-mode projections or by any generated plane
state.

## Bounds and fixed design

Each coordinate halfwidth is one half of the observed score range among the 18
validated forcing-generated B landing states:

```text
|alpha_1| <= 0.0406483232,
|alpha_2| <= 0.0705250541.
```

No positivity shrink was required. Across the full rectangle, the minimum
saved concentrations are `min(u) = 0.330184` and `min(v) = 2.48394`.

The frozen 17-by-17 master grid contains 289 complete states:

- 25 points in a common 5-by-5 initial seed;
- 16 points in an interleaved, untouched 4-by-4 fixed test set;
- 248 remaining acquisition-pool points.

Every manifest row has `simulation_status = not_run` and an empty outcome
label. Thus this build cannot leak a B/C label into later model selection.

## Reflection audit

For every state, the output retains the raw, reflected, feature-canonical, and
serpentine-continuity-oriented fields. Raw and reflected inputs canonicalize to
the same state to numerical precision.

The feature-canonical path contains a genuine orientation discontinuity: its
largest step is about `94.3` times the corresponding continuity-oriented step.
Consequently, later plots must retain a continuity-oriented diagnostic and must
not interpret the feature-canonical jump as a fold or gap in the physical
sampling family.

## Files and verification

The construction is implemented by:

- `brusselator_ml_boundary_geometry_config.m`;
- `build_ml_BC_boundary_geometry.m`;
- `validate_ml_BC_boundary_geometry.m`.

The output is in
`experiment_outputs/ml_BC_boundary_geometry_v1/`. It contains the complete
states, manifest, source provenance, bracket-projection audit, bound sources,
reflection-continuity audit, figure, and machine-readable configuration.

All 16 frozen-geometry validation checks pass. The next computation is to evolve
the 25 initial-seed states with the frozen PDE and classify them as B, C, other,
or unresolved using the exact orbit library and adaptive finite-time horizon.
The exact edge orbit remains reserved for evaluating the learned uncertainty
set and dynamically recovered brackets.
