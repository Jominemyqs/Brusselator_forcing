# Dynamic B--C edge tracking: recurrent boundary-state candidate

## Purpose

The multiclass slice work supplies a finite-time B--C initial-condition bracket

```text
lambda_B = 0.187500000,
lambda_C = 0.187890625.
```

Both endpoints lie on the same reflection-aligned ramp-40/ramp-45 landing-state
line and differ by `2.11e-4` in the symmetric physical full-state norm. A long
trajectory near a static bisection point is not by itself an edge trajectory,
so `run_dynamic_edge_tracking_BC` implements an actual evolve--separate--
rebracket cycle.

## Dynamic algorithm

The two exact frozen `b = 10` trajectories are evolved in their maintained
common spatial orientation. Rebracketing is triggered when their symmetric
quadrature-weighted full-state separation reaches `2e-2`. At the common event
time, the complete evolved `(u,v)` states are linearly interpolated and every
midpoint is classified against A/B/C/R5/R6.

Classification is adaptive: trajectories are first tested for 80 time units
and extended in 40-unit chunks, up to 500, without changing the established
late median or maximum thresholds. The bracket is updated only for a unique B
or C label. Any A, R5, R6, ambiguous, or unresolved midpoint stops the search.

The initial pair reaches the separation trigger after 66.6 time units. Seven
midpoint classifications reduce the evolved-state bracket to

```text
alpha_B = 0.8046875,
alpha_C = 0.8125000,
state separation = 1.6493e-4.
```

No third outcome occurs. The midpoint outcomes and first sustained commitment
times are

| Step | Alpha | Outcome | Decision time |
|---:|---:|---:|---:|
| 1 | 0.5000000 | B | 133 |
| 2 | 0.7500000 | B | 165 |
| 3 | 0.8750000 | C | 193 |
| 4 | 0.8125000 | C | 237 |
| 5 | 0.7812500 | B | 181 |
| 6 | 0.7968750 | B | 200 |
| 7 | 0.8046875 | B | 248 |

The renewed pair remains below the same separation trigger for another 92
time units. One rebracketing therefore increases the assembled boundary-
shadowing interval from 66.6 to 158.6 time units.

The arithmetic midpoint of close B/C states is saved only as an approximate
visual shadow. It is not treated as an exact PDE trajectory. All recurrence
claims below are computed separately from the two exact bracketing
trajectories.

## Recurrent dynamics

`analyze_dynamic_edge_BC_candidate` analyzes the exact 92-unit
post-rebracketing segment. Both sides repeatedly cross the same mean-`v`
Poincare section and exhibit a common period near 4.003:

| Window start | Side | Period | Period CV | One-cycle median return | One-cycle maximum return |
|---:|:---:|---:|---:|---:|---:|
| 0 | B | 4.00297 | `4.31e-5` | `1.34e-4` | `2.99e-4` |
| 0 | C | 4.00370 | `4.19e-4` | `3.16e-4` | `1.42e-3` |
| 20 | B | 4.00294 | `4.62e-5` | `1.13e-4` | `1.67e-4` |
| 20 | C | 4.00415 | `4.32e-4` | `3.91e-4` | `1.42e-3` |
| 40 | B | 4.00286 | `4.66e-5` | `9.62e-5` | `1.38e-4` |
| 40 | C | 4.00492 | `4.29e-4` | `5.61e-4` | `1.42e-3` |

The candidate remains separated from the validated orbit library. Across the
two exact trajectories, its smallest observed phase/reflection-aware distance
to any known outcome is `0.1954`; the minima to B are about `0.2209`, and the
minima to C are `0.1954`--`0.1990`. Its period also differs from B
(`3.9483`), C (`4.3113`), R5 (`5.0585`), and R6 (`6.5837`).

## Transverse growth

Sampling the B/C separation once per B-side Poincare return removes the strong
within-cycle oscillation. For separations below `1e-2`, 17 samples fit

```text
d(t) proportional to exp(gamma t),
gamma = 0.0538981,
R^2 = 0.9999979.
```

The corresponding finite-time separation multiplier over one candidate period
is `exp(gamma T) = 1.24079`. Fits truncated at separations `5e-3`, `1e-2`, and
`1.5e-2` give effectively the same rate and multiplier.

This exponential growth is consistent with a transverse unstable direction.
It is not a Floquet multiplier calculation: the trajectory has not yet been
converged to an exact periodic solution, and the measured B/C difference need
not be a pure eigenvector.

## Multicycle convergence and fixed phase

`run_dynamic_edge_tracking_BC_multicycle` continues the original event with
three additional evolve--separate--rebracket cycles. The three renewed
segments last `88.2`, `88.2`, and `91.8` time units before reaching the same
`2e-2` separation threshold. All 21 new midpoint classifications are B or C,
and the assembled edge-shadow time reaches `426.8`.

The original per-segment mean-v section is useful for recurrence but allows
its centering constant to drift. `analyze_dynamic_edge_BC_fixed_phase` therefore
rephases all four segments on the single Newton-compatible hyperplane

```text
<X-X_ref,F(X_ref)>_w = 0
```

with positive crossing derivative. It phase-fixes the exact B and C trajectories
separately before averaging their section states. Every side has one regular
upward branch per period; the minimum normalized transversality is approximately
one. In the final event the period is `4.00282885`, the B/C phase-state distance
is `9.57e-5`, and the direct one-period residual of the averaged seed is
`2.63e-4`. Consecutive fixed-phase template differences remain at a numerical
floor of roughly `2.6e-4`--`5.5e-4`, so the analysis exports a shooting seed
rather than claiming convergence from edge tracking alone.

## Newton-converged periodic orbit

`run_periodic_edge_orbit_newton` solves

```text
Phi_T(X)-X = 0,
<X-X_ref,F(X_ref)>_w = 0
```

by weighted, matrix-free finite-difference Newton--GMRES shooting. There is no
spatial-shift variable because homogeneous Neumann boundaries leave reflection
as a discrete symmetry but remove continuous spatial translation. Three
accepted Newton corrections give

| Quantity | Value |
|---|---:|
| Period | `4.0029149136` |
| Tight-ode15s normalized shooting residual | `1.17e-11` |
| Fixed phase-condition residual | `1.61e-18` |
| Independent tight-ode45 flow residual | `2.15e-9` |
| State correction from edge seed | `2.38e-4` |

This verifies a periodic orbit of the stated `N = 400`, `Lx = 40`, frozen
`b = 10` semidiscrete system.

## Floquet spectrum

`run_periodic_edge_orbit_floquet` integrates the exact variational equation
along the converged orbit and applies the monodromy map through weighted
matrix-free Arnoldi iteration. The leading results are

| Mode | Multiplier | Modulus | Real Floquet exponent |
|---|---:|---:|---:|
| Nontrivial unstable | `1.24226303` | `1.24226303` | `0.0541942` |
| Temporal phase | `1.00000000` | `1.00000000` | approximately zero |
| Leading stable real mode | `0.64374` | `0.64374` | `-0.11003` |
| Leading stable complex pair | `0.39459 +/- 0.40242 i` | `0.56360` | `-0.14325` |

All ten reported eigenpair residuals are below `2.5e-9`. The tangent action for
the unstable eigenvector agrees with a nonlinear finite-difference flow-map
check to relative error `2.13e-6`. Exactly one computed nontrivial multiplier
lies outside the unit circle. Its difference from the earlier finite-time
separation estimate `1.24079` is only `0.001473`.

The saved `periodic_edge_orbit_floquet_v2` result supersedes `v1`: the first run
computed the real phase and unstable modes correctly but reapplied the real
monodromy routine directly to complex eigenvectors when auditing complex-mode
residuals. Version 2 applies the real operator separately to real and imaginary
parts; the multipliers themselves are unchanged and all residual audits pass.

## Unstable branches

`run_periodic_edge_orbit_unstable_branches` perturbs the exact orbit on both
sides of the real unstable Floquet eigenvector. The established five-way
phase/reflection-aware classifier gives

| Relative amplitude | Negative branch | Positive branch |
|---:|:---:|:---:|
| `1e-5` | B | C |
| `1e-4` | B | C |

Decision times are `240` and `273` at amplitude `1e-5`, and `197` and `231` at
amplitude `1e-4`, respectively. No branch is forced into a binary label: A,
R5, R6, ambiguity, and unresolved behavior remain allowed by the classifier.

## Current conclusion

For the `N = 400`, `Lx = 40` semidiscrete frozen system, the evidence chain is

```text
dynamic B/C rebracketing
  -> phase-audited recurrent shooting seed
  -> Newton-converged periodic orbit
  -> exactly one nontrivial unstable Floquet multiplier
  -> opposite unstable branches approach B and C.
```

It is now supported to write:

> The local B--C basin boundary is organized by an unstable periodic edge
> state with one unstable Floquet direction.

This statement is local and initially discretization-qualified. The subsequent
fixed-domain exact-orbit continuation at `N = 200, 400, 600` retains the orbit
and exactly one nontrivial unstable multiplier at every grid. Aligned full-orbit
differences decrease by a factor `5.48` from the coarse--middle comparison to
the middle--fine comparison, consistent with the second-order stencil. See
`PERIODIC_EDGE_RESOLUTION_RESULTS.md`.

That remaining mechanism gate is now complete. B and C are separately
Newton-converged and Floquet-attracting at `N = 200, 400, 600`, and both signs
of each consistently oriented edge eigenvector reach opposite exact B/C orbit
libraries. Thus the complete local B--edge--C geometry, not only the edge orbit,
persists across the tested grids. See `BC_EDGE_RESOLUTION_RESULTS.md`.

The narrow C--B--A slice interval near `lambda = 0.276` remains deliberately
unmapped at high resolution. It is reserved as a held-out discovery test for
the proposed uncertainty-guided ML method; these results still do not justify
direct C--A edge tracking or a claim of fractal basin geometry.

## Reproducible outputs

- Driver: `run_dynamic_edge_tracking_BC.m`
- Configuration: `brusselator_dynamic_edge_BC_config.m`
- Adaptive classifier: `brusselator_evolve_and_classify_state.m`
- Candidate analysis: `analyze_dynamic_edge_BC_candidate.m`
- Multicycle driver: `run_dynamic_edge_tracking_BC_multicycle.m`
- Fixed-section analysis: `analyze_dynamic_edge_BC_fixed_phase.m`
- Newton driver: `run_periodic_edge_orbit_newton.m`
- Floquet driver: `run_periodic_edge_orbit_floquet.m`
- Unstable-branch driver: `run_periodic_edge_orbit_unstable_branches.m`
- Numerical output: `experiment_outputs/dynamic_edge_tracking_BC_v1`
- Analysis output: `experiment_outputs/dynamic_edge_BC_candidate_analysis_v1`
- Multicycle output: `experiment_outputs/dynamic_edge_tracking_BC_multicycle_v1`
- Fixed-phase output: `experiment_outputs/dynamic_edge_BC_fixed_phase_convergence_v1`
- Periodic-orbit output: `experiment_outputs/periodic_edge_orbit_newton_v1`
- Floquet output: `experiment_outputs/periodic_edge_orbit_floquet_v2`
- Branch output: `experiment_outputs/periodic_edge_orbit_unstable_branches_v1`
- Resolution driver: `run_periodic_edge_orbit_resolution_continuation.m`
- Resolution output: `experiment_outputs/periodic_edge_orbit_resolution_continuation_v1`
- Matched B--edge--C driver: `run_BC_edge_resolution_geometry.m`
- Matched-geometry output: `experiment_outputs/BC_edge_resolution_geometry_v1`
- Borderline-branch confirmation: `confirm_BC_edge_resolution_borderline_branch.m`
