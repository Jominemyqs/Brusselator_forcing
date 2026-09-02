# Dynamic edge tracking from the held-out C--B--A region

## Scope

This experiment tests the two exceptionally tight state-space brackets obtained
by the post-budget, five-way refinement of the frozen held-out challenge:

```text
C--B: lambda in [0.2761230469, 0.2761840820]
B--A: lambda in [0.2761840820, 0.2762451172]
width of each bracket = 6.103515625e-5
```

The original eight-run held-out acquisition did not sample B. These brackets
therefore remain post-budget follow-up objects, not successful frozen-budget ML
discoveries. Their dynamical tests are nevertheless useful: the C--B side has
an independently computed reference edge orbit, while the B--A side tests
whether endpoint labels continue to define a direct two-basin boundary.

All evolving-state rebracketing uses the unchanged five-way physical classifier
(`A`, `B`, `C`, `R5`, `R6`, with `unresolved` retained), a physical pair-distance
restart target of `2.5e-4`, and a separation event at pair distance `0.02`.
Unexpected outcomes stop a nominally pairwise run rather than being coerced to
one endpoint.

## C--B result: recovery of the known periodic edge state

The first dynamic cycle narrowed the original C--B pair with eight midpoint
labels and produced a restart separation of `1.33798e-4`. That pair separated
again after only `1.0` time unit and remained far from the known edge orbit,
with minimum orbit distance about `0.22366`. This first failure shows why a
single rebracketing event is not enough to identify the organizing object.

Three further dynamic rebracketing cycles then gave:

| Segment | Initial separation | Time to separation `0.02` | Minimum distance to exact `E_BC` | Time inside `d(E_BC)<0.02` |
|---:|---:|---:|---:|---:|
| 2 | `1.9415e-4` | `89.8` | `3.9082e-4` | `61.4` |
| 3 | `1.5706e-4` | `88.2` | `3.5546e-5` | `88.2` |
| 4 | `1.5768e-4` | `92.0` | `9.6437e-5` | `92.0` |

All 21 new midpoint integrations resolved to C or B; no third outcome or final
unresolved label obstructed the sequence. The nearly constant later shadow
durations and the `3.55e-5` minimum physical orbit distance show that repeated
dynamic rebracketing recovers the independently Newton-converged periodic edge
orbit `E_BC` (previously measured period `4.0029149` and leading nontrivial
Floquet multiplier `1.242263`). In segments 3 and 4, the tracked exact side
trajectories remain within the documented `0.02` edge-orbit neighborhood for
the complete retained segment.

The supported conclusion is:

> The C--B boundary intersecting the narrow held-out C--B--A state-space region
> is dynamically connected to the stable manifold of the already verified
> periodic edge state `E_BC`.

This is the desired dynamical-usefulness test: a boundary bracket localized in
the challenge region leads, through PDE edge tracking, to an invariant object
computed independently of the ML classifier.

## Nominal B--A result: five-way obstruction by C

The first B--A cycle completed eight B/A midpoint labels, reached restart
separation `1.71565e-4`, and shadowed for `3.2` time units. At the next evolved
state-space chord, however, the guarded refinement produced:

| Step | Local interpolation coordinate | Outcome | Tested time |
|---:|---:|---:|---:|
| 1 | `0.500` | B | `120` |
| 2 | `0.750` | B | `80` |
| 3 | `0.875` | C | `160` |

The run therefore stopped with status `event_obstructed_by_C`. The evolved B/A
endpoint pair does not define a simple direct B--A transition along its current
linear chord. This is not a failed numerical bisection; it is a multiclass
geometric result and validates the decision to retain the five-way classifier
during dynamic edge tracking.

No B--A edge state or direct B--A basin boundary is established by this run.
The appropriate continuation is a moderately dense five-way scan of the
evolved B--A chord, followed by edge tracking only on genuinely adjacent label
intervals. In particular, the new C intrusion must not be silently assigned to
either B or A.

## Reproducibility

Implementation:

- `brusselator_ml_CBA_edge_tracking_config.m`
- `run_ml_CBA_dynamic_edge_tracking.m`
- `brusselator_ml_CBA_multicycle_config.m`
- `run_ml_CBA_dynamic_edge_multicycle.m`

Machine-readable outputs and figures:

- `experiment_outputs/ml_heldout_CBA_edge_tracking_CB_v1`
- `experiment_outputs/ml_heldout_CBA_edge_tracking_CB_multicycle_v1`
- `experiment_outputs/ml_heldout_CBA_edge_tracking_BA_v1`
- `experiment_outputs/ml_heldout_CBA_edge_tracking_BA_multicycle_v1`

The multicycle result files store every midpoint classification, the evolved
endpoint trajectories, restart separations, segment durations, exact-edge
distances where applicable, configuration, and run metadata. Output directories
are versioned and the runners refuse to overwrite completed result files.
