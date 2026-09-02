# Five-way scan after the nominal B--A edge obstruction

## Question

After one dynamic B--A rebracketing cycle, a midpoint on the evolved endpoint
chord reached C. The endpoint labels therefore no longer justified binary B--A
bisection. This scan asks what outcome ordering is actually sampled along that
evolved chord.

## Design

The exact B and A endpoint states are loaded from the first completed dynamic
cycle at global edge time `7.6`. The scanned family is

```text
X(alpha) = (1-alpha) X_B + alpha X_A,  0 <= alpha <= 1.
```

The nonuniform grid uses spacing `0.1` through `alpha=0.7` and spacing `0.025`
from `0.75` through `1`. The five previously determined states at `alpha = 0,
0.5, 0.75, 0.875, 1` are reused exactly. Fourteen new states are evolved under
the frozen `b=10` PDE and classified against A/B/C/R5/R6 with the unchanged
physical phase/reflection-aware classifier and an adaptive horizon through at
most `t=500`.

The endpoints at `alpha=0` and `1` retain their B and A labels by forward
invariance of their already established basin outcomes. They are not counted
as new PDE labels, and no late-distance values are attached to those propagated
endpoint labels.

## Result

All 19 sampled states receive known outcomes. The sampled ordering is

```text
B: alpha = 0, 0.1, ..., 0.7, 0.75, 0.775, 0.8, 0.825
C: alpha = 0.85, 0.875
A: alpha = 0.9, 0.925, 0.95, 0.975, 1
```

There are no R5, R6, or final unresolved outcomes. The two candidate transition
intervals at this resolution are

| Left outcome | Right outcome | Interval | Width |
|---|---|---:|---:|
| B | C | `[0.825, 0.850]` | `0.025` |
| C | A | `[0.875, 0.900]` | `0.025` |

The C label at `alpha=0.85` is not marginal: it resolves after a long transient
with decision time `168`, tested duration `200`, and late median physical
C-orbit distance `4.79e-3`. The A cases at and above `0.9` resolve rapidly,
with decision times `13`--`18`. The 14 new integrations record about `1037.3 s`
of solver runtime, dominated by the C case at `alpha=0.85`.

The supported conclusion is deliberately resolution-qualified:

> At the sampled resolution, the dynamically evolved chord that was nominally
> B--A crosses a finite C-basin region and has the ordering B--C--A. It cannot be
> interpreted as a direct binary B--A slice.

Neighboring different labels remain candidate adjacent intervals, not proof
that no still narrower outcome window lies between them.

## Bounded next step

The B--C organizing mechanism is already independently established as the
periodic edge orbit `E_BC`, including exact Newton convergence, Floquet spectrum,
opposite B/C unstable branches, and grid continuation. Repeating that complete
pipeline is not the current bottleneck.

The only scientifically necessary local closure is the new C--A side:

1. refine `[0.875,0.9]` with the five-way guard to approximately `1e-3`;
2. stop refinement immediately if B, R5, R6, or a persistent unresolved state
   appears;
3. if the bracket remains clean C/A, run a bounded dynamic rebracketing pilot
   of at most three rebracketing events under the current 220-unit segment cap,
   to determine whether it approaches reproducible recurrent boundary dynamics
   or is obstructed after evolution;
4. do not begin Newton/Floquet continuation of a new object unless that pilot
   first supplies a reproducible candidate.

That is the recommended stopping milestone for the current project phase.

## Closure outcome

The bounded closure has now been run. Its first midpoint at `alpha=0.8875`
reaches B by tested frozen time `200`, with late median B-orbit distance
`2.5544e-3`. The prespecified guard therefore stops the experiment before any
dynamic C--A rebracketing. At the newly sampled resolution the local ordering
is B--C--B--A, and no direct C--A edge is established. See
`BOUNDED_CA_CLOSURE_RESULTS.md`.

## Reproducibility

- Configuration: `brusselator_ml_evolved_BA_scan_config.m`
- Driver: `run_ml_evolved_BA_five_way_scan.m`
- Output: `experiment_outputs/ml_heldout_CBA_evolved_BA_scan_v1`

The output contains all new raw trajectories, reused-state provenance, the
complete label table, transition and region tables, figures, metadata, and a
restart checkpoint. The driver refuses to overwrite a completed result.
