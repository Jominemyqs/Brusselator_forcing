# First full-state edge-tracking attempt: the B--A line intersects C

## Question

The first edge pilot asked whether the physically realized ramp-40 B landing
and ramp-45 A landing provide a direct two-basin bracket under the frozen
`b=10` dynamics.

The complete landing fields were reflection-aligned and interpolated as

`X(lambda) = (1-lambda) X_B + lambda X_A`,

with `lambda=0` the B landing and `lambda=1` the A landing. Every trajectory
was classified in the quadrature-weighted physical full-state norm against
the validated A/B/C/R5/R6 library. A midpoint was allowed to update a binary
bracket only after a unique A or B classification.

## Result

The fresh endpoint integrations reproduce B at `lambda=0` and A at
`lambda=1`. The first midpoint, `lambda=0.5`, reaches A with a sustained
decision time near 7. The second midpoint, `lambda=0.25`, does not reach A or
B. It approaches C and first satisfies the sustained C-neighborhood test near
time 119.

At the end of the moderate 240-unit run, the C distance is `1.19e-4`; its
late median and maximum are `3.26e-4` and `3.11e-3`. The other late median
distances remain between approximately 0.34 and 0.55. This is convergence
toward the known C orbit, not evidence for an edge trajectory.

The result was repeated from exactly the same initial state with tight
`ode15s` tolerances (`RelTol=1e-8`, `AbsTol=1e-10`, `MaxStep=0.1`). The tight
run again reaches C, with late median `2.76e-4`, late maximum `3.05e-3`, and
final distance `9.71e-5`. The moderate and tight C-distance histories nearly
coincide.

## Interpretation

The 40/45 forcing protocols have different final outcomes, but their landing
states do not provide a direct B--A basin-boundary bracket along this linear
state-space segment. The C basin intervenes:

- `lambda=0`: B;
- `lambda=0.25`: C;
- `lambda=0.5`: A;
- `lambda=1`: A.

Therefore an A--B binary bisection on this line would be scientifically
incorrect. The point at `lambda=0.25` is also not a C-accessing forcing
protocol: it is a constructed frozen-system initial condition. It establishes
that C is accessible from the state-space geometry between two physical
landings, while C remains unreached by the sampled forcing family itself.

## Revised next step

The line must be treated as a multibasin slice. The next computations should
refine two pairwise transitions independently:

1. the B--C transition inside `lambda in (0,0.25)`;
2. the C--A transition inside `lambda in (0.25,0.5)`.

Each refinement must retain the five-way stopping rule because another basin
could intervene. Only after a genuinely adjacent pairwise bracket is found
should repeated dynamic rebracketing be used to shadow an edge trajectory.
The two boundaries may have different organizing edge states.

This result strengthens the multistable interpretation: even a line joining
physical B and A landings can pass through the basin of a frozen attractor that
the forcing scan has not directly reached.

## Follow-up status

The follow-up scan and guarded refinements are reported in
`EDGE_SLICE_MULTICLASS_RESULTS.md`. They establish a longer finite-time B--C
bracket, but invalidate the apparent C--A adjacency by finding a narrow B
outcome at finer resolution. The independent C--A bisection proposed above is
therefore superseded by a local multiclass remap.
