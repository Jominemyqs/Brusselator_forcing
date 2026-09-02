# Project status and recommended temporary conclusion

## Current research question

The project is now best organized around:

> How can finite-rate parameter forcing produce an apparently binary protocol
> threshold while the nearby forcing-return region intersects several frozen
> basins, which invariant structures organize those multibasin intersections,
> and can task-aware probabilistic sampling locate dynamically useful brackets
> for classical edge tracking?

This is narrower and better supported than claiming a universal overshoot
threshold, a direct one-to-one correspondence between protocol thresholds and
frozen basin boundaries, or generic ML discovery of invariant objects.

## What is established numerically

For the documented `Lx=40` problem:

1. The frozen `b=10` system has a stationary patterned outcome A and several
   distinct periodic patterned outcomes. B and C have exact numerical periodic
   representations and local Floquet evidence; R5 and R6 have grid/solver,
   recurrence, cross-preparation, and sampled local-attraction evidence.
2. The tested physical forcing family reaches A, B, R5, and R6-type recurrent
   regimes. It has not reached C. This is sampled non-reachability, not proof
   that C is inaccessible.
3. The original threshold near `bmax=11.13` is reproduced in its historical
   setting but is seed- and resolution-sensitive. It should not be the main
   universal claim.
4. The local B--C separator is mechanistically resolved: a Newton-converged
   periodic edge orbit has one nontrivial unstable Floquet direction, and its
   two unstable branches reach B and C. The matched B--edge--C geometry persists
   at `N=200,400,600` for fixed `Lx=40`.
5. One-dimensional state-space slices are genuinely multiclass. The original
   A--B slice reaches C, the held-out challenge contains a narrow C--B--A
   structure, and the newly evolved nominal B--A chord samples B--C--A.
6. The first pointwise uncertainty acquisition improves proper scores and
   separator-band uncertainty but fails to propose a useful B/C bracket. A
   pair-aware acquisition succeeds retrospectively in the controlled B--C
   plane, recovering `E_BC` from three of four proposed pairs. In the frozen
   held-out challenge, the eight-run endpoint misses the narrow B basin; a
   separately reported post-budget refinement finds it and again recovers
   `E_BC`. General rare-basin discovery is therefore not yet established.
7. A bounded, five-way refinement of the physical symmetric-ramp family at
   `bmax=11.14`, `Thold=80`, and seed 1 used 20 protocol runs. It found only B
   and A among the sampled protocols and localized their finite-time change to
   `T_ramp in (40.6005859375, 40.6015625)`, an interval of width
   `9.765625e-4`. This is a sampled, configuration-specific threshold, not a
   universal critical rate. No sampled protocol reached C, R5, or R6.
8. The two adjacent B/A forcing-return states are already far apart in the
   physical full-state norm (`d approximately 0.535`), and each is already in
   its final attractor neighborhood when `b` returns to 10. Thus they do not
   supply a close frozen B/A edge bracket: most of the selection has occurred
   during the nonautonomous forcing phase.
9. Guarded interpolation of those adjacent return states is nevertheless
   multibasin. A C midpoint interrupted B/A rebracketing and was independently
   confirmed with the tighter solver. The locally adjacent B/C pair recovered
   the exact `E_BC` orbit under dynamic edge tracking, with minimum orbit
   distance `4.15e-4` and 61 time units inside the predefined exact-edge
   neighborhood.

## Temporary scientific conclusion

A defensible current conclusion is:

> Finite-rate overshoots act as structured probes of a multistable patterned
> state space. A one-parameter forcing family can display an extremely sharp,
> apparently binary B/A selection threshold even though a chord through its
> neighboring return states intersects the C basin. The B/C part of that local
> multibasin geometry is organized by the robust periodic edge state `E_BC`.
> Hence protocol adjacency, state-space adjacency, and basin adjacency are
> distinct notions. Pair-aware probabilistic proposals can recover dynamically
> meaningful edge brackets in controlled settings, but direct discovery of
> rare forcing-accessible basins and the invariant object governing the
> physical B/A protocol threshold remain open.

This conclusion has both a dynamical-systems contribution and an honest ML
result. It does not identify `E_BC` as the cause of the physical B/A threshold:
C was not reached by the sampled forcing protocols, and the adjacent return
states were already strongly separated at the end of forcing. It also does not
claim that every observed outcome is globally robust, that C is physically
unreachable, or that ML learned an invariant manifold.

## Recommended stopping milestone

Complete one bounded closure experiment on the newly exposed C--A candidate:

1. five-way refine the sampled interval `[0.875,0.9]` to roughly `1e-3`;
2. if it remains C/A, perform at most three dynamic rebracketing events under
   the current 220-unit segment cap;
3. stop when the pilot either produces reproducible recurrent boundary
   dynamics or encounters a third outcome/unresolved obstruction.

Then pause computation and write the results around the temporary conclusion
above. Do not proceed automatically to a new Newton--Krylov solve, Floquet
analysis, full four-parameter forcing atlas, autoencoder, Koopman model, or
full-state active learner. Each is a separate next research phase and should be
justified by what the bounded C--A pilot reveals.

## Stopping milestone reached

The first guarded midpoint at `alpha=0.8875` reaches B rather than C or A. The
closure therefore stops after one new trajectory, exactly as prespecified, and
no dynamic C--A pilot is run. The sampled evolved chord now has ordering
B--C--B--A. This reinforces the temporary conclusion that thin or folded
multiclass basin intersections—not a single binary separator—organize this
state-space region. It does not establish fractal, riddled, or Wada geometry.

The current computational phase is complete. The next activity should be
manuscript organization and figure selection. Further refinement of the new
C--B or B--A subintervals belongs to a separately scoped basin-topology phase,
not to the present closure experiment.

## Physical forcing closure reached

The additional physical forcing experiment has now also reached its bounded
stopping condition:

1. all 11 prespecified half-unit seed protocols on `[40,45]` were run;
2. the remaining nine-run budget refined the only sampled label change;
3. the final physical interval is B/A with width `9.765625e-4`;
4. no forcing-generated C outcome was sampled;
5. a third outcome C stopped binary frozen rebracketing of the return states;
6. the C result passed a tighter-solver confirmation;
7. the resulting B/C chord pair recovered the independently known `E_BC`.

This phase should stop here. A denser search for a very narrow forcing-accessible
C window, computation of a nonautonomous B/A edge object, or continuation of a
new B/A separator would each be a separately declared research phase rather
than an unreported extension of the closure budget.
