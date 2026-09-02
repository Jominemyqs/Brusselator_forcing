# Exact periodic edge-orbit resolution continuation

## Question

Does the Newton-converged B--C periodic edge state and its one-dimensional
instability persist when the spatial grid is changed at fixed

```text
L = 40, b = 10, a = 2.5, d2 = 9.73, sigma = 0.45?
```

This study changes only the method-of-lines resolution. It is not a domain-size
study or a proof of a periodic solution of the continuum PDE.

## Method

`run_periodic_edge_orbit_resolution_continuation` reuses the exact `N = 400`
orbit, regrids its phase-fixed state to `N = 200` and `N = 600` with
shape-preserving cubic interpolation, and recomputes the target-grid frozen
vector field for the phase condition. Each target then undergoes the same
tight single-shooting Newton--GMRES solve and independent tight-`ode45`
verification as the source orbit.

The leading Floquet spectrum is recomputed from the exact target-grid
variational equation. Complete periodic trajectories are resampled on a common
spatial/temporal grid and compared after minimizing over temporal phase and the
exact spatial reflection.

## Exact-orbit and Floquet results

| N | dx | Period | Newton residual | Independent `ode45` residual | Unstable multiplier | Real exponent |
|---:|---:|---:|---:|---:|---:|---:|
| 200 | `0.201005` | `4.0017303566` | `5.50e-10` | `2.60e-9` | `1.2495914933` | `0.0556801` |
| 400 | `0.100251` | `4.0029149136` | `1.17e-11` | `2.15e-9` | `1.2422630263` | `0.0541942` |
| 600 | `0.066778` | `4.0031465171` | `4.70e-12` | `2.21e-9` | `1.2408623676` | `0.0539092` |

At every resolution:

- Newton shooting converges;
- the independent integrator reproduces the orbit to below `3e-9`;
- the phase multiplier differs from one by at most `1.32e-9`;
- exactly one computed nontrivial multiplier lies outside the unit circle;
- every reported leading eigenpair residual is below `2.8e-9`;
- the tangent/nonlinear finite-difference check is approximately `2.2e-6`.

## Full-orbit convergence

| Grid pair | Period difference | Phase-fixed state distance | Aligned full-orbit distance | Selected phase shift | Reflection |
|:---:|---:|---:|---:|---:|:---:|
| 200--400 | `1.1846e-3` | `5.6548e-3` | `6.0968e-3` | `0` | no |
| 400--600 | `2.3160e-4` | `1.0345e-3` | `1.1124e-3` | `0` | no |

The full-orbit difference ratio is `5.4806`. This is close to the ratio expected
when differences are proportional to differences in `dx^2` on these three
unequally spaced grids. The comparison also shows that the common phase
condition remains consistent: neither a temporal shift nor reflection is needed
to align the continued orbits.

A diagnostic least-squares fit of the form `Q(dx) = Q_0 + c dx^2` gives

```text
period extrapolate       4.0033161482
unstable-multiplier extrapolate 1.2398023563
```

The maximum three-point fit residuals are `6.52e-6` for the period and
`2.47e-5` for the multiplier. These extrapolates support second-order numerical
consistency but should not be presented as rigorous continuum limits from only
three grids.

## Conclusion and next gate

The exact periodic edge orbit and its single nontrivial unstable direction
persist across `N = 200, 400, 600` at fixed `L = 40`. This substantially reduces
the risk that the mechanism found at `N = 400` is a grid-specific artifact.

That next mechanistic robustness test is complete. B and C are exact and
Floquet-attracting on the same grids, and the consistently oriented unstable
edge branches select B/C at every resolution. See
`BC_EDGE_RESOLUTION_RESULTS.md`.

## Reproducible files

- Configuration: `brusselator_periodic_edge_resolution_config.m`
- Driver: `run_periodic_edge_orbit_resolution_continuation.m`
- Reusable Newton driver: `run_periodic_edge_orbit_newton.m`
- Reusable Floquet driver: `run_periodic_edge_orbit_floquet.m`
- Output: `experiment_outputs/periodic_edge_orbit_resolution_continuation_v1`
