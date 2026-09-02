# Adaptive frozen-PDE labeling of the initial B--C seed

## Purpose and classifier

`ml_BC_boundary_initial_labeling_v1` evolves the 25 predeclared initial-seed
states from `ml_BC_boundary_geometry_v1` under the actual frozen `b = 10` PDE
at `N = 400` and `Lx = 40`.

The five-way physical classifier uses:

- the independently Newton-converged exact N=400 B and C periodic orbits;
- the validated Stage-2 A, R5, and R6 references;
- the unchanged quadrature-weighted full-state distance, minimized over exact
  reflection and temporal phase for periodic references;
- late-median and late-maximum thresholds `0.01` and `0.02`.

Reported labels are mapped without nearest-class coercion:

```text
B -> B
C -> C
A, R5, R6 -> other
unresolved or ambiguous -> U
```

`U` means unresolved at the maximum finite horizon; it is not an attractor.
The exact edge orbit and held-out C--B--A challenge region are not used by the
labeling library.

## Adaptive time horizon

Every state is first integrated to `t = 80`. A trajectory that has not passed
a unique physical outcome threshold is extended in 40-unit chunks, up to
`t = 400`. Per-case partial files and a study-level checkpoint make this
procedure resumable without changing any label definitions.

No case resolved at `t = 80`. The final tested-duration distribution was:

| Tested duration | Number of states |
|---:|---:|
| 120 | 5 |
| 160 | 13 |
| 200 | 5 |
| 240 | 1 |
| 280 | 1 |

This confirms that adaptive extension is necessary even for apparently clear
points in this initial plane.

## Outcomes

All 25 states eventually received unique B or C labels:

```text
B:     14
C:     11
other:  0
U:      0
```

Arranged by increasing `alpha_2` from bottom to top and increasing `alpha_1`
from left to right, the 5-by-5 label map is

```text
B B B B B
B B B B B
C C B B B
C C C C B
C C C C C
```

Thus the coarse separator intersection is curved. In particular,
`(alpha_1,alpha_2)=(0.0406483,0.0352625)` reaches B while the corresponding
`alpha_2` level reaches C in the preceding four columns. The grid contains
seven horizontal or vertical nearest-neighbor pairs with opposite B/C labels.

The plane midpoint took until `t = 280` to satisfy the B threshold, with first
sustained commitment at approximately `t = 200`. The upper-right point took
until `t = 240` to satisfy the C threshold, committing at approximately
`t = 224`. These long decision times provide an independent finite-time signal
of difficult near-separator regions, although decision time alone is not a
basin-boundary definition.

## Reproducibility and limits

The complete output occupies approximately 141 MB because every full trajectory
is retained. It includes 25 per-case raw files, the exact label mapping, outcome
library provenance, progress and final tables, configuration metadata, and a
figure reproducible from the saved result using
`plot_ml_BC_boundary_initial_labels.m`.

All 14 post-run validation checks pass. These results establish labels only on
the coarse 25-point seed and do not yet reconstruct the separator between those
points. The next step is to fit the first calibrated probabilistic B/C model,
use uncertainty to rank the frozen acquisition pool, and compare its proposed
queries with a budget-matched uniform strategy. The 16-point fixed test set
remains unlabeled and untouched.
