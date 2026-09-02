# Matched B--edge--C resolution geometry

## Scientific question

Does the complete local basin-boundary mechanism

```text
attracting B orbit <- unstable periodic edge orbit -> attracting C orbit
```

persist when the spatial resolution changes at fixed `L = 40` and `b = 10`?

The preceding study continued only the edge orbit. This study separately
Newton-converges B and C, verifies their leading Floquet spectra, constructs
exact grid-matched B/C reference libraries, and repeats both unstable edge
branches at `N = 200, 400, 600`.

## Exact attracting periodic orbits

| Outcome | N | Period | Newton residual | Largest nontrivial Floquet modulus |
|:---:|---:|---:|---:|---:|
| B | 200 | `3.9505240544` | `6.45e-12` | `0.874605` |
| B | 400 | `3.9482993142` | `7.21e-12` | `0.875885` |
| B | 600 | `3.9478937452` | `1.03e-11` | `0.876115` |
| C | 200 | `4.3211863638` | `8.74e-12` | `0.840047` |
| C | 400 | `4.3113478763` | `4.43e-11` | `0.843415` |
| C | 600 | `4.3095272291` | `2.75e-10` | `0.844057` |

All independent tight-`ode45` residuals are below `2.7e-9`. No computed
nontrivial multiplier lies outside the unit circle, and all leading eigenpair
residuals are below `3.9e-9`. This supports local attraction of both exact
periodic orbits throughout the tested grid range.

The aligned complete-orbit differences are

| Outcome | Pair | Full-orbit distance | Selected phase shift | Reflection |
|:---:|:---:|---:|---:|:---:|
| B | 200--400 | `6.6609e-4` | `0` | no |
| B | 400--600 | `1.2593e-4` | `0` | no |
| C | 200--400 | `5.4262e-3` | `0` | no |
| C | 400--600 | `1.0109e-3` | `0` | no |

The corresponding decrease factors are `5.29` for B and `5.37` for C,
consistent with the second-order spatial discretization on these unequal grids.

## Consistent unstable direction

Floquet eigenvectors have arbitrary signs, so each target-grid unstable edge
direction is oriented by its physical weighted inner product with the regridded
`N = 400` direction. The resulting overlaps are

```text
N=200: 0.99915
N=400: 1
N=600: 0.99997
```

Thus the unstable direction itself persists smoothly under grid refinement.

## Grid-matched branch outcomes

Each exact edge orbit is perturbed at relative physical amplitude `1e-4` and
integrated under the frozen system. Classification minimizes physical distance
over temporal phase and exact reflection using the exact target-grid B and C
orbits. The classifier contains an unresolved fallback; it does not choose the
nearest orbit when neither tolerance is met.

| N | Negative branch | Positive branch | Latest decision time |
|---:|:---:|:---:|---:|
| 200 | B | C | `225` |
| 400 | B | C | `231` |
| 600 | B | C | `231` |

At `N = 200`, the positive branch first meets the C criterion at global time
240 with late median/maximal C distances `9.31e-3` and `1.68e-2`. Because this
is close to the declared thresholds, an additional fixed 80-unit continuation
was performed. At global time 320 its late median/maximal exact-C distances are
`6.43e-4` and `2.29e-3`, confirming continued approach rather than a brief
neighborhood crossing.

## Conclusion

For fixed `L = 40`, the matched local geometry is robust across the tested
spatial resolutions:

```text
exact attracting B
    <- negative unstable branch
exact periodic edge state with one unstable multiplier
    -> positive unstable branch
exact attracting C.
```

This is stronger than persistence of the edge orbit alone. It supplies a
dynamically verified, resolution-tested separator against which a learned
uncertainty boundary can be evaluated.

## Reproducible files

- Configuration: `brusselator_BC_edge_resolution_config.m`
- Main driver: `run_BC_edge_resolution_geometry.m`
- Borderline confirmation: `confirm_BC_edge_resolution_borderline_branch.m`
- Main output: `experiment_outputs/BC_edge_resolution_geometry_v1`
- Confirmation output: `experiment_outputs/BC_edge_resolution_branch_confirmation_v1`
