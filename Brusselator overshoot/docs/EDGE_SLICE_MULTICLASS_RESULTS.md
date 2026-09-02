# Multiclass landing-state slice scan and guarded refinement

## Question

The first ramp-40/ramp-45 full-state interpolation reached B at
`lambda = 0`, C at `lambda = 0.25`, and A at `lambda = 0.5`. Endpoint labels
therefore did not justify independent binary searches without first checking
for additional outcome intervals. This study maps the constructed slice with
the unchanged physical A/B/C/R5/R6 classifier before dynamic edge tracking.

These are frozen-system initial-condition experiments at `b = 10`. Except for
the two physical endpoints, intermediate states are linear combinations of
forcing landing states and are not themselves produced by a forcing protocol.

## Numerical definitions

- `Lx = 40`, `N = 400`;
- tight-moderate `ode15s`, output interval and `MaxStep` 0.2,
  `RelTol = 1e-6`, `AbsTol = 1e-8`;
- physical quadrature-weighted full-state distances minimized over temporal
  phase and exact reflection;
- late median threshold `1e-2` and late maximum threshold `2e-2`;
- initial frozen duration 80, extended to 240 only when unresolved;
- every midpoint is classified against A, B, C, R5, and R6. No binary bracket
  is updated for a third, ambiguous, or unresolved outcome.

## Moderately dense five-way scan

`run_five_way_edge_slice_scan` evaluates 21 points at `Delta lambda = 0.025`
over `0 <= lambda <= 0.5`. It exactly reuses the saved `lambda = 0`, `0.25`,
and `0.5` pilot trajectories. The sampled regions are

| Sampled interval | Outcome | Samples |
|---|---:|---:|
| `0 <= lambda <= 0.175` | B | 8 |
| `0.20 <= lambda <= 0.275` | C | 4 |
| `0.30 <= lambda <= 0.50` | A | 9 |

No R5, R6, ambiguous, or unresolved label remains at the maximum tested time.
At this resolution the only candidate transition intervals are
`(0.175, 0.20)` for B--C and `(0.275, 0.30)` for C--A. They are candidate
adjacencies at the sampled resolution, not proofs that no narrower outcome
region intervenes. Outputs are in `experiment_outputs/edge_slice_scan_BCA_v1`.

## Multiclass-guarded refinements

`run_edge_slice_pairwise_refinement` applies six guarded midpoint tests to
each candidate interval.

### B--C side

| `lambda` | Outcome | First sustained decision time |
|---:|---:|---:|
| 0.187500000 | B | 181 |
| 0.193750000 | C | 157 |
| 0.190625000 | C | 173 |
| 0.189062500 | C | 189 |
| 0.188281250 | C | 208 |
| 0.187890625 | unresolved at 240 | -- |

The increasing decision times show slow near-boundary dynamics, but do not
identify an invariant edge state. At `t = 240`, the last midpoint has a late
median C distance `1.65e-2` and a finite-window recurrence candidate near
period `4.2433`. Its recurrence alone is not a new recurrent class because the
trajectory is drifting toward the validated C orbit.

`run_edge_slice_BC_unresolved_extension` continues exactly that trajectory to
`t = 500`. It commits to C at about `t = 241`; its final-window median,
maximum, and final C distances are `2.35e-4`, `2.17e-3`, and `3.36e-4`.
Accordingly, the longer finite-time B--C bracket is

```text
0.187500000 (B) < lambda* < 0.187890625 (C),
width = 0.000390625.
```

This is an initial-condition outcome bracket, not a computed stable manifold
or edge state.

### Putative C--A side

The refinement first gives rapid A outcomes down to `lambda = 0.2765625`, and
C at `lambda = 0.27578125`. However, their midpoint
`lambda = 0.276171875` reaches B at about `t = 86`. Therefore the apparent
C--A adjacency is invalidated at finer resolution. At minimum, the local
sampled order is

```text
C at 0.275781250, B at 0.276171875, A at 0.276562500.
```

This may be a narrow B tongue, an interleaved basin structure, or only the
first resolved part of a still finer multiclass arrangement. The present data
do not justify direct C--A edge tracking.

Outputs are in `experiment_outputs/edge_slice_pairwise_refinement_v1` and
`experiment_outputs/edge_slice_BC_unresolved_extension_v1`.

## Scientific interpretation

The constructed landing-state line intersects all three validated basins, but
not as a single monotone B--C--A partition at all tested scales. A moderately
dense scan sees three broad regions; guarded refinement then exposes a much
narrower B outcome inside the apparent C--A transition. This is precisely why
endpoint-only binary bisection is unsafe in the multistable PDE.

The B--C side exhibits long but ultimately C-directed transients. The current
evidence supports a sharp finite-time outcome boundary and supplies a close
B/C initial-condition pair, but it does not yet show a reproducible invariant
object on that boundary. The C--A side must be locally remapped before any
edge-state claim or pairwise bisection.

The result does not change the forcing-accessibility statement: the C outcomes
on this slice are constructed initial conditions, while the sampled physical
forcing protocols still reach A, B, R5, or R6 but not C.

## Next computations

1. Run a local five-way scan across
   `[0.27578125, 0.2765625]`, retaining the C, B, and A samples already
   computed, to determine whether clean C--B and B--A adjacent intervals exist.
2. Refine only adjacencies supported by that local scan, again with five-way
   stopping rules.
3. Develop dynamic edge tracking first for the longer-time B--C bracket
   `[0.1875, 0.187890625]`: periodically rebracket full evolving states rather
   than interpreting a single long transient as the edge trajectory.
4. Compare any sustained boundary dynamics with the forcing landing curve
   only after numerical persistence and recurrence have been established.
