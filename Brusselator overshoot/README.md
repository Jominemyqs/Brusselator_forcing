# 1D Brusselator overshoot experiments

This repository contains the MATLAB simulation code, Python probabilistic-
acquisition code, frozen experiment definitions, and compact scientific audit
records supporting the manuscript
*Forcing Accessibility, Periodic Edge States, and Task-Aware Probabilistic
Acquisition in the One-Dimensional Brusselator*.

> **Release status.** This is the first version-controlled research snapshot.
> Raw numerical trajectories remain outside ordinary Git history because the
> local output tree is several gigabytes. Before public release, the selected
> data package must be deposited in a persistent archive and its DOI recorded
> in `CITATION.cff` and `docs/REPRODUCIBILITY.md`.

## Quick navigation

- [`docs/REPOSITORY_GUIDE.md`](docs/REPOSITORY_GUIDE.md): logical code map and
  principal scientific workflows;
- [`docs/REPRODUCIBILITY.md`](docs/REPRODUCIBILITY.md): environment, data, and
  verification instructions;
- [`docs/MANUSCRIPT_DATA_MANIFEST.csv`](docs/MANUSCRIPT_DATA_MANIFEST.csv):
  checksums for the compact result tables used by the manuscript Supplement;
- [`docs/CNSNS_RELEASE_CHECKLIST.md`](docs/CNSNS_RELEASE_CHECKLIST.md): code and
  data tasks to complete before journal submission;
- [`docs/README.md`](docs/README.md): index of experiment-design and result
  notes.

## Quick start

Python analyses use Python 3.9 or newer:

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements-dev.txt
python -m pytest
python scripts/verify_repository.py
```

MATLAB experiments are run from the repository root so their declared
`experiment_outputs/...` paths remain valid. A no-simulation configuration
check is:

```matlab
cfg = brusselator_default_config();
brusselator_validate_config(cfg);
```

The legacy course-project scripts and their saved results are retained unchanged.
New reproducible studies use the small configuration layer introduced in July 2026.

## Repository layout

- MATLAB and Python files at the repository root implement simulations,
  analysis, plotting, and acquisition benchmarks.
- `docs/` contains versioned experiment designs, audits, and result notes.
- `experiment_outputs/` contains the current versioned numerical studies.
- Legacy top-level `*_outputs/` directories and `*_results.mat` files are kept
  in their historical locations so old scripts remain reproducible.
- The manuscript, current paper figures, literature library, and archived
  course-paper material live one level above this repository in the surrounding
  `Overshoot/` workspace.

Large generated data are intentionally excluded from Git by `.gitignore` while
remaining available locally. Scientific conclusions and compact result records
belong in `docs/`; use `git add -f` only when a generated artifact has been
deliberately selected for version control.

## Running a new study

From this directory in MATLAB, run one of:

```matlab
run_robustness_study
run_high_resolution_stiff_check
run_solver_crosscheck
run_extended_postforcing_trajectories
analyze_phase_aware_attractors
run_frozen_full_state_confirmation
analyze_full_state_poincare_recurrence
run_periodic_orbit_local_stability
analyze_local_stability_orbit_distance
run_extended_even_sector_stability
run_phase_aware_basin_slice
extend_phase_aware_basin_slice_edge_candidates
analyze_edge_candidate_dynamics
run_edge_candidate_frozen_confirmation
run_candidate_c_local_stability
run_candidate_c_numerical_validation
run_forcing_to_c_reachability
run_ramp_transition_refinement
run_unresolved_forcing_long_extensions
analyze_unresolved_recurrent_groups
analyze_recurrent_group_threshold_sensitivity
run_recurrent_class_validation
run_recurrent_class_local_stability
analyze_recurrent_class_local_stability_dense
analyze_three_way_phase_aware_remap
run_dynamic_edge_tracking_BC_multicycle
analyze_dynamic_edge_BC_fixed_phase
run_periodic_edge_orbit_newton
run_periodic_edge_orbit_floquet
run_periodic_edge_orbit_unstable_branches
run_periodic_edge_orbit_resolution_continuation
run_BC_edge_resolution_geometry
confirm_BC_edge_resolution_borderline_branch
```

Each script creates a new subdirectory of `experiment_outputs/` and refuses to
overwrite an existing study directory. It saves the configuration, MATLAB/runtime
metadata, summary CSV/MAT data, and raw numerical data separately.

`brusselator_default_config.m` is the single source for the model, grid, forcing,
initial-condition, solver, and recovery defaults. Its final integration time is
always computed as `Tup + Thold + Tdown + post_forcing`, preventing the prior
`Thold = 320, Tfinal = 300` scheduling error.

The robustness screen uses `N = [200, 400, 800]`, seeds `[1, 2]`, MATLAB-default
and tight `ode45` tolerances, and the bracketing amplitudes 11.12 and 11.14. Its
separate `run_high_resolution_stiff_check` script provides a tight-`ode15s`
cross-check at `N = 800`; it has no maximum-step constraint and therefore serves
as a preliminary solver-resolution indicator, rather than a convergence result.
`run_solver_crosscheck` makes the corresponding `N = 400` comparison with
`MaxStep = 0.1`. The extended-trajectory study uses tight `ode15s` with a
documented maximum step and integrates the 11.12 and 11.14 cases for 10,000 time
units after forcing ends. Run `plot_extended_postforcing_results(outdir)` to
recreate its figures from saved raw trajectories only.
These are finite-time tests; they do not establish an asymptotic attractor by
themselves.

The phase-aware analysis reads the long-trajectory output only. It minimizes the
late-time `v`-field distance over temporal lag and reflection, and records
profile/spectral statistics. The frozen full-state confirmation restarts the
saved baseline and departure states at `b = 10`, with fine temporal output,
and saves both fields. The multi-return analysis then tests direct and
reflection-related Poincare returns from those saved full states.

`run_phase_aware_basin_slice` is a deliberately different replacement for the
legacy final-`v` slice scripts. It interpolates the complete `(u,v)` state between
the late stationary frozen-control state `A` and a Poincare-standardized state
`B` on the candidate periodic orbit. Frozen `b = 10` trajectories are compared
to `A` modulo reflection and to the entire orbit modulo both temporal phase and
reflection. It records an unresolved outcome rather than forcing a binary basin
label. `extend_phase_aware_basin_slice_edge_candidates` continues only the
unresolved saved cases to a longer final time with the identical metric.
`analyze_edge_candidate_dynamics` then performs a late-window full-state
Poincare/recurrence and spectral analysis of the persistent unresolved case;
it is a characterization step, not an attractor label.

## Current generated studies

`experiment_outputs/robustness_reference_v2` contains 18 completed, checkpointed
members of the 24-case `ode45` matrix. The remaining high-resolution `ode45`
members are intentionally left pending because they are costly; the two separate
stiff-solver cross-checks are saved alongside it. The saved results show that the
historical threshold is reproduced at `N = 400`, seed 1, but changes with both
seed and resolution. It must therefore be regarded as a non-robust finite-time
observation until a controlled convergence study resolves that dependence.

The current `N = 400`, seed-1 frozen-state evidence identifies a stationary
baseline pattern and a distinct candidate time-periodic patterned state after
the `bmax = 11.14` overshoot. In a 200-unit full-state restart at `b = 10`
with output spacing and maximum step `0.05`, the departure trajectory returns
close to its spatial reflection after one section crossing (about `1.974` time
units) and close to itself after two crossings. This is finite-time evidence
for a spatiotemporally reflection-symmetric periodic state. A local stability
screen then perturbed three standardized phases in two smooth Neumann-compatible
symmetry sectors at relative amplitudes `1e-4` and `1e-3`; all 12 trajectories
returned to the periodic-orbit neighborhood. The initially slow reflection-even
response was extended to 160 time units (about 40 direct periods): its `1e-5`
and `1e-4` perturbations returned to the unperturbed restart floor, and its
`1e-3` perturbation decayed toward the same scale. This supports local attraction in
the sampled directions, but is not a Floquet calculation, a full basin proof,
or a robust continuum result.

The first phase-aware slice (`phase_aware_basin_slice_v1`) confirms the endpoint
distinction in this one `N = 400`, seed-1 frozen-system setting: `lambda <= 0.725`
returned to the stationary baseline by time 80, whereas the standardized periodic
endpoint `lambda = 1` returned to the periodic-orbit neighborhood. The intermediate
`lambda = 0.75` and `0.8` cases were unresolved at time 80. The targeted extension
to time 240 found that `lambda = 0.75` approaches the periodic orbit (late median
orbit distance about `3.2e-4`), but `lambda = 0.8` remains far from both endpoint
neighborhoods (late medians about `0.56` to the baseline and `0.45` to the orbit).
At that point, this was finite-time evidence for slow or more complicated
transitional dynamics along the slice, not evidence for a third attractor or a
computed basin boundary.

The late-window diagnostic for the `lambda = 0.8` trajectory changes that last
interpretation in a specific but still limited way. Over frozen times 160--240,
it has 18 upward spatial-mean-`v` Poincare crossings with a crossing-period
median about `4.3099` and coefficient of variation `3.4e-4`. Consecutive full
Poincare states return directly with median relative distance `2.7e-4` (and are
not reflection-related); longer return distances drift gradually but remain below
`5e-3` through 16 crossings. Its late dominant mean-field frequency is about
`0.2372` (period about `4.22` at the saved-output resolution). These observations
are finite-time evidence for a distinct candidate periodic patterned state `C`.
The fresh frozen restart in `edge_candidate_frozen_confirmation_v1` gives the
same direct recurrence in a separate 160-unit restart segment: its last 80 units have
period `4.3113`, period CV `6.0e-5`, and consecutive-return distance `1.82e-4`
(with a best nine-crossing distance `2.47e-5`). The matched 12-case local screen
then returned to C in every tested phase, parity sector, and amplitude
(`1e-4` or `1e-3`) over 40 frozen time units; late median orbit distances were at
most `1.18e-3`, and late maxima at most `1.04e-2`. Together these provide
finite-time support that C is locally attracting in the sampled directions.
They still do not prove a third stable attractor, its full basin, or robustness
to grid, seed, domain, and solver changes.

`run_candidate_c_numerical_validation` is the next gate. It transfers a late
candidate-C state to several grids, a second time integrator, a looser but still
controlled tolerance set, and modestly changed domains. A case supports C only
when it has close full-state Poincare recurrence, regular section timing, a
period close to the source orbit, and a phase-standardized spatial profile close
to C. The domain transfer is performed in normalized `x/L`, so those cases test
persistence after deformation rather than independent discovery. Independent
preparation is reserved for the subsequent forcing-to-C reachability scan.
Domain sensitivity is reported separately from the fixed-`Lx` numerical gate:
changing `Lx` changes the PDE problem and may cross a patterned-state
bifurcation, whereas changing `N`, solver, or tolerances tests the numerical
representation of the same `Lx = 40` problem.

`run_forcing_to_c_reachability` runs only after the fixed-`Lx` validation gate.
It saves both fields for a focused amplitude/duration/rate/asymmetry protocol
sample, classifies late frozen trajectories as A/B/C/unresolved, adaptively
extends unresolved cases, and repeats any apparent C hit with the tight solver.

## Candidate-C validation and forcing reachability (August 2026)

`candidate_c_numerical_validation_v1` supports C as a numerically credible
finite-time periodic outcome for the stated `Lx = 40` problem. Regridded
restarts at `N = 200`, `600`, and `800` all retained a direct period near
`4.31`, close full-state recurrence, and the standardized C profile. The
`N = 200` `ode45` cross-check agreed with `ode15s`, and the moderate-tolerance
`N = 400` result agreed with the tight reference. The modest domain transfers
did not retain the same C identification: `Lx = 38` moved toward a roughly
`1.99`-period recurrent regime and `Lx = 42` toward a distinct roughly
`4.68`-period regime. Claims about C must therefore remain domain-qualified.

`forcing_to_c_reachability_pilot_v1` tested 20 full-state protocols at
`N = 400`, `Lx = 40`: an amplitude line, hold durations from 0 to 320,
symmetric ramp times from 10 to 160, and two ramp asymmetries. It found five A
outcomes, twelve B outcomes, no C outcomes, and three unresolved trajectories
after the adaptive horizon (the `bmax = 12.5`, `bmax = 13`, and symmetric
10-unit-ramp cases). Thus C is independently validated as a frozen-system
outcome but has not been prepared by any tested physical overshoot protocol.
This is a finite sampled non-reachability result, not proof that C is
unreachable. The most informative next scan is a refinement of the ramp-time
transitions between unresolved/B (`10 < T_ramp < 20`) and B/A
(`20 < T_ramp < 80`), using landing-state geometry rather than another broad
parameter grid.

`ramp_transition_refinement_v1` resolves the two transitions on the symmetric
ramp-time line at `bmax = 11.14` and `Thold = 80` without changing the
three-way A/B/C definitions. Ramp times 10--16 remain unresolved after 240
post-forcing time units, 17--40 approach B, and 45--80 approach A. The sampled
transition brackets are therefore `(16,17)` between the unresolved regime and
B, and `(40,45)` between B and A. No tested ramp reaches C.

`unresolved_forcing_long_extension_v1` continues every unresolved forcing case
to a common post-forcing time of 1000. None of the nine cases approaches A, B,
or C over that window. Instead, all nine exhibit regular Poincare timing and
close full-state recurrence. Two cases (`bmax = 12.5` and the 16-unit symmetric
ramp) have period about `5.0585`; the other seven have period about `6.5837`.
The best late recurrence distances are between about `3.4e-5` and `1.5e-4`.
These are reproducible finite-time recurrent-pattern candidates rather than
unstructured unresolved transients, but they are not yet established stable
attractors.

`unresolved_recurrent_group_analysis_v1` compares one-period full-state
templates after minimizing over temporal phase and spatial reflection. A
single exploratory cutoff of `0.02` produces five candidate groups, including
two groups reached by multiple protocols. The companion
`unresolved_recurrent_threshold_sensitivity_v1` shows why five must not be
reported as an attractor count: the grouping changes from nine to two connected
components as the cutoff increases from `0.005` to `0.04`. The robust statement
is the existence of two strongly separated recurrence-period classes. Their
minimum cross-class full-state distance is about `0.738`, whereas the maximum
within-class distances are about `0.0246` and `0.0758`. Whether the spatial
variants inside the `6.5837` class are distinct invariant orbits remains open.

Together these studies complete the first finite-time ramp-resolution stage.
They also show that the A/B/C classifier is incomplete for the forcing-generated
landing set. Subsequent landing-state geometry should retain the original
A/B/C labels while displaying these recurrent cases as candidate additional
outcomes; it should not silently fold them into C or treat them as one proven
new attractor.

## Recurrent-class validation and Stage 2 design (August 2026)

`recurrent_class_validation_v1` applies a reduced version of the candidate-C
gate to representatives of the two recurrence-period classes. The temporary
names are R5 (period about `5.0585`, prepared by `bmax = 12.5`) and R6 (period
about `6.5837`, prepared by the 12-unit symmetric ramp). Periods are confirmed
by independent, phase-fixed, full-state frozen restarts rather than inferred
only from Poincare crossing times. All eight grid/solver cases support periodic
orbit persistence: tight `ode15s` at `N = 200, 400, 600` and an `ode45`
cross-check at `N = 200` for each class. Direct return errors are approximately
`1.0e-3`--`2.4e-3`, and the section-period CVs are below `3e-5`.

The tight cross-preparation test sends both R5 preparations to the same orbit
set. Four of seven R6 preparations meet the fixed `0.02` same-orbit cutoff by
the end of the 80-unit confirmation. The remaining ramp-13, ramp-14, and
ramp-15 cases retain the R6 period and move substantially toward its reference,
but finish at distances `0.0201`, `0.0250`, and `0.0227`; they therefore remain
provisional R6 preparations rather than being silently relabeled. Both classes
remain about `0.75` apart in the full-orbit metric.

The initial 360-phase local-stability analysis exposed an amplitude-independent
nearest-phase floor in its maximum-distance statistic. The versioned
`recurrent_class_local_stability_dense_v1` reanalyzes the same saved PDE
trajectories against 1440-phase `pchip` templates. All four tested local
perturbations per class return to their representative orbit neighborhoods.
For R5, late median distances are `5.4e-4`--`6.6e-4` and late maxima are below
`9.91e-3`; for R6, late medians are about `4.5e-4` and late maxima are below
`6.53e-3`. The sampled directions use reflection-odd and reflection-even smooth
Neumann-compatible perturbations at relative amplitudes `1e-4` and `1e-3`.

Taken together, these tests provide finite-time numerical support that R5 and
R6 are distinct locally attracting periodic outcomes for the stated `Lx = 40`
problem. This wording remains intentionally narrower than a Floquet stability
calculation, a full basin result, a continuum convergence theorem, or a domain
robustness claim.

The next projection stage is specified in `docs/STAGE2_LANDING_GEOMETRY_DESIGN.md`
and centralized by `brusselator_stage2_geometry_config.m`. It uses two views:
a phase-fixed representation of A/B/C/R5/R6 and forcing landing states, and a
full-orbit representation in which periodic outcomes appear as loops. The
primary POD basis is balanced across the five validated outcomes; landing
states are projected passively, and provisional/novel/unresolved labels remain
available. The 40/45 B--A and 16/17 R5--B brackets are carried explicitly into
the state manifest for subsequent edge tracking.

`three_way_phase_aware_remap_v1` reclassifies saved states against all three
references. Using full `(u,v)` late-window distances for the frozen slice,
the initial 80-unit runs give A for `lambda <= 0.725`, B for `lambda = 1`,
and unresolved labels for `lambda = 0.75` and `0.8`; their 240-unit extensions
give B at `lambda = 0.75` and C at `lambda = 0.8`. This demonstrates that the
straight A--B slice is a three-outcome diagnostic, not a binary basin test.
The historical 10,000-unit overshoot output saved only `v(x,t)`, so its
late-window labels are explicitly V-only and agree with its saved full final
snapshots: control and `bmax = 11.12` are A, while `bmax = 11.14` is B. None of
these three saved overshoot trajectories approaches C. A future forcing scan
must save both fields throughout to support fully phase-aware trajectory labels.

## Stage 2 landing-state geometry (August 2026)

`build_stage2_landing_geometry` implements the specification in
`docs/STAGE2_LANDING_GEOMETRY_DESIGN.md`. The finalized versioned output is
`stage2_landing_geometry_v2`; the first diagnostic build remains preserved as
`stage2_landing_geometry_v1` and is not the reported result. The v2 manifest
contains 39 unique forcing protocols and verifies that every saved landing is
the complete finite `(x,u,v)` state at the exact time that `b` returns to 10 on
the declared `N = 400`, `Lx = 40` grid.

The formal phase-section audit passes for B, C, R5, and R6. B has two upward
mean-`v` crossings per period, but the two branches differ by only
`2.03e-4` in the physical full-state norm after exact reflection reduction, so
they are treated as reflection-equivalent. C, R5, and R6 each have one upward
crossing. Crossing counts persist under temporal and `N = 200` spatial
resampling, all normalized transverse derivatives exceed 1.16, and the late
cycle intersections reproduce the corresponding templates to below
`8.22e-4`.

The primary POD uses one stationary A sample and 80 phase samples from each
periodic outcome, with explicit aggregate weight `0.2` for each of A, B, C,
R5, and R6. No copies of A are stored. The first two, three, and five modes
explain `69.47%`, `78.45%`, and `90.92%` of this balanced outcome covariance.
Removing the separate `u/v` RMS normalization leaves the landing geometry
qualitatively stable: pairwise landing-distance correlations are 0.996, 0.992,
and 0.994 in two, three, and five dimensions. These are visualization
diagnostics only; all reported orbit distances continue to use the physical
quadrature-weighted full-state norm in the original variables.

Reflection handling is materially relevant. Fourteen of 39 landing states use
the reflected feature-canonical orientation, and several ordered forcing
families show larger artificial steps after independent feature
canonicalization than under the saved nearest-neighbor continuity orientation.
The raw, reflected, feature-canonical, and family-continuity-oriented states
are all retained. Canonicalizing a raw state or its reflected copy agrees to a
maximum relative error of `6.08e-9`.

The physical landing distances sharpen the two edge-tracking targets. The
16-unit symmetric-ramp landing is assigned to R5 and is `0.0141` from the R5
orbit, whereas the 17-unit landing is assigned to B and is `0.00945` from B.
The 40-unit landing is `0.00401` from B, while the 45-unit landing is `0.00214`
from A. No forcing landing is close to C: the smallest sampled C distance is
`0.377`, from the 20-unit hold protocol. This remains sampled evidence that C
is not accessed by the present forcing family, not an inaccessibility proof.

The orbit-only POD has a 20-mode median landing reconstruction error of
`0.0270`, but its worst error is `0.511`, at the 16-unit R5--B endpoint. That
exceeds the largest orbit-library error by a factor of 8.00 and triggers the
predeclared secondary landing-aware POD diagnostic. With half the aggregate
weight retained on the balanced orbit library and half distributed over the
39 landings, its 20-mode median and maximum landing errors fall to `0.00543`
and `0.0637`. This secondary basis does not replace the primary coordinates or
define a basin boundary; it shows that linear landing-aware coordinates are
adequate before introducing nonlinear representation learning.

Stage 2 therefore supports moving next to full-state edge tracking on the
physically realized 40/45 B--A bracket, followed by the 16/17 R5--B bracket.
The POD plots organize initial guesses and expose symmetry artifacts, but the
edge computation must continue to classify trajectories with physical
phase/reflection-aware orbit distances rather than POD proximity.

## First B--A edge-tracking attempt (August 2026)

`edge_tracking_BA_landing_v1` applies a five-outcome-safe full-state bisection
to the reflection-aligned ramp-40 B and ramp-45 A landing states. Fresh frozen
integrations confirm the endpoints, and the midpoint at `lambda = 0.5` reaches
A. The next midpoint, `lambda = 0.25`, reaches C rather than A or B: its late
median C distance is `3.26e-4`, late maximum is `3.11e-3`, and final distance
is `1.19e-4` after 240 time units. It first remains inside the prescribed C
neighborhood near time 119.

`edge_tracking_BA_landing_C_confirmation_v1` repeats exactly that constructed
initial state with tight `ode15s` tolerances. It again reaches C, with late
median `2.76e-4`, late maximum `3.05e-3`, and final distance `9.71e-5`. The
moderate and tight distance histories nearly coincide.

Consequently, the 40/45 physical forcing outcomes do not define a direct B--A
basin bracket along their linear landing-state interpolation: the C basin
intervenes. The `lambda = 0.25` state is constructed and is not a new
C-accessing forcing protocol. The next edge work must refine the B--C interval
inside `(0,0.25)` and the C--A interval inside `(0.25,0.5)` separately, retaining
the five-way stopping rule. See `docs/EDGE_TRACKING_BA_RESULTS.md` for the method,
evidence, and revised interpretation.

## Multiclass slice scan and guarded refinement (August 2026)

`edge_slice_scan_BCA_v1` performs the required five-way scan before binary
edge bisection: 21 reflection-aligned full-state initial conditions at
`Delta lambda = 0.025` over `0 <= lambda <= 0.5`, with adaptive extension from
80 to 240 frozen time units. The sampled sequence is B on `[0,0.175]`, C on
`[0.20,0.275]`, and A on `[0.30,0.50]`. No sampled R5, R6, ambiguous, or
persistently unresolved case remains. These intervals describe the sampled
line only and do not establish basin adjacency.

`edge_slice_pairwise_refinement_v1` then refines both candidate transitions
while retaining the five-way stopping rule. On the B--C side, decision times
grow to roughly 181--208 time units as the interval narrows. The midpoint at
`lambda = 0.187890625` is still unresolved at 240, but the separate
`edge_slice_BC_unresolved_extension_v1` continuation commits to C near time
241 and is close to C through time 500. The resulting longer finite-time
B--C bracket is `(0.1875,0.187890625)`, of width `3.90625e-4`.

The apparent C--A interval is not binary. The refinement finds C at
`lambda = 0.27578125`, A at `0.2765625`, and B at their midpoint
`0.276171875`. Thus a narrow B region, or a still finer interleaved structure,
was invisible to the `0.025` scan. Direct C--A edge tracking is not justified;
the local interval must first be remapped with a finer multiclass scan. The
B--C long transient is also not itself an edge orbit: repeated dynamic
rebracketing is the next requirement. See `docs/EDGE_SLICE_MULTICLASS_RESULTS.md`
for the full finite-time evidence and revised handoff.

## Dynamic B--C edge-state candidate (August 2026)

`dynamic_edge_tracking_BC_v1` performs the first genuine dynamic
evolve--separate--rebracket cycle on the longer-time B--C bracket. The initial
pair begins `2.11e-4` apart and reaches the declared `2e-2` separation trigger
after 66.6 frozen time units. Seven multiclass-safe midpoint classifications
restore the evolved pair to separation `1.65e-4`; no A, R5, R6, ambiguous, or
unresolved obstruction appears. The renewed pair shadows common dynamics for
another 92 time units before reaching the trigger again, giving 158.6 assembled
shadowing time units.

`dynamic_edge_BC_candidate_analysis_v1` finds the same recurrent dynamics on
both exact sides of the renewed bracket. The mean period is approximately
`4.00354`; period CVs are below `4.32e-4`, and one-cycle median return errors
are `9.62e-5`--`5.61e-4` across the tested windows. The trajectories remain at
least `0.195` from every validated A/B/C/R5/R6 outcome.

At B-side Poincare returns, the B/C separation follows an exponential fit with
growth rate `0.053898`, `R^2 = 0.999998`, and an implied finite-time multiplier
`1.24079` per candidate period. These results support an unstable periodic
B--C edge-state candidate near period 4.0035.

The subsequent `dynamic_edge_tracking_BC_multicycle_v1` run adds three dynamic
rebracketing events without encountering A, R5, R6, or another obstruction.
It extends the assembled shadow to `426.8` frozen time units. Reanalysis on one
fixed Newton-compatible phase hyperplane gives one transverse upward crossing
per cycle on both exact sides, relative transversality near one, a final period
`4.00282885`, and a best seed flow residual `2.63e-4`.

`periodic_edge_orbit_newton_v1` turns that seed into a Newton--GMRES-converged
single-shooting solution at

```text
T = 4.0029149136
ode15s normalized flow residual = 1.17e-11
independent tight-ode45 residual = 2.15e-9
```

`periodic_edge_orbit_floquet_v2` then finds exactly one leading nontrivial
unstable multiplier, `mu_u = 1.24226303`, the trivial temporal-phase multiplier
at one, and all other computed leading multipliers inside the unit circle. The
unstable value is within `0.001473` of the independent edge-separation estimate
`1.24079`. Finally, `periodic_edge_orbit_unstable_branches_v1` perturbs the orbit
along both signs of its unstable eigenvector: at relative amplitudes `1e-5` and
`1e-4`, the negative branch reaches B and the positive branch reaches C.

For the stated `N = 400`, `Lx = 40`, frozen-`b = 10` semidiscrete problem, this
supports the mechanistic conclusion that the local B--C basin boundary is
organized by a periodic edge state with one unstable Floquet direction. Grid,
solver, and domain continuation of the exact orbit and spectrum remain required
before interpreting it as a continuum-robust PDE object. See
`docs/DYNAMIC_EDGE_BC_RESULTS.md`.

`periodic_edge_orbit_resolution_continuation_v1` completes the first fixed-domain
grid continuation. Regridded `N = 200` and `N = 600` seeds both Newton-converge,
and all three resolutions retain exactly one nontrivial unstable multiplier:

| N | Period | Unstable multiplier |
|---:|---:|---:|
| 200 | `4.00173036` | `1.24959149` |
| 400 | `4.00291491` | `1.24226303` |
| 600 | `4.00314652` | `1.24086237` |

The phase/reflection-aligned full-orbit difference decreases from `6.10e-3`
between `N = 200` and `400` to `1.11e-3` between `N = 400` and `600`. No
temporal shift or reflection is selected in either comparison. The difference
ratio `5.48` and the nearly linear period/multiplier trends in `dx^2` are
consistent with the second-order finite-difference stencil. The three-grid
extrapolates (`T approximately 4.00331615`, `mu_u approximately 1.23980236`)
are diagnostics, not continuum proofs. See `docs/PERIODIC_EDGE_RESOLUTION_RESULTS.md`.

`BC_edge_resolution_geometry_v1` completes the matched local-geometry test.
B and C are Newton-converged separately at `N = 200, 400, 600`; both have no
computed nontrivial unstable multiplier at any grid. Their largest nontrivial
Floquet moduli remain near `0.876` for B and `0.844` for C. The edge unstable
eigenvectors have physical overlaps `0.99915`, `1`, and `0.99997` with the
oriented `N = 400` direction. At relative amplitude `1e-4`, every matched-grid
test gives

```text
-v_u -> B,
+v_u -> C.
```

The exact B/C libraries contain only those two outcomes but retain an unresolved
fallback, so a nonmatching trajectory is never forced into a binary label. The
initially marginal `N = 200` positive branch was extended to global time 320;
its late median exact-C distance decreases from `9.31e-3` to `6.43e-4`. See
`docs/BC_EDGE_RESOLUTION_RESULTS.md`.

This verified separator is now a suitable controlled ground truth for the first
ML experiment. The proposed experiment uses only two landing-aware POD
coordinates, a probabilistic B/C classifier with explicit other/unresolved
handling, uncertainty-guided PDE sampling, and edge-orbit recovery as the main
metric. The narrow `lambda approximately 0.276` C--B--A interval is held out as
a discovery test. See `docs/ML_BOUNDARY_DISCOVERY_DESIGN.md`.

`ml_BC_boundary_geometry_v1` now freezes that experiment's pre-simulation
geometry. The reflection-aligned midpoint of the clean B/C bracket is used as
the center, and the first two landing-aware POD modes are decoded into physical
`u,v` directions. A 17-by-17 master grid provides 25 initial-seed points, 16
untouched test points, and 248 acquisition candidates. All 289 states pass the
positivity and coordinate-reprojection gates, and none has an inherited or
simulated outcome label. The exact edge orbit and the narrow C--B--A challenge
region are excluded from construction. The two-mode plane retains `97.59%` of
the known bracket displacement, while the reflection audit exposes a real
canonical-orientation jump that must remain a plotting diagnostic. See
`docs/ML_BC_BOUNDARY_GEOMETRY_RESULTS.md`.

`ml_BC_boundary_initial_labeling_v1` then evolves the frozen 25-point initial
seed using exact N=400 B/C periodic references and validated A/R5/R6 references.
The unchanged physical classifier is evaluated first at time 80 and unresolved
cases are extended in 40-unit chunks to at most time 400. The final seed has 14
B and 11 C outcomes, with no `other` or unresolved cases. No trajectory resolved
at the first 80-unit horizon; the longest required 280 units. The resulting
5-by-5 map contains seven nearest-neighbor B/C brackets and visibly curves near
positive `alpha_1`. The exact edge and held-out C--B--A interval remain excluded,
and the fixed test set remains unlabeled. See
`docs/ML_BC_BOUNDARY_LABELING_RESULTS.md`.

`ml_BC_probabilistic_boundary_v1` fits the first boundary model using only those
25 labels. It uses a deterministic ARD squared-exponential probit GP with
Laplace inference; nested leave-one-out refitting gives 96% accuracy, and the
single error is the curved-boundary B point at positive `alpha_1`. The inferred
`p_C = 0.5` contour is a single curved segment. Eight uncertainty/novelty-guided
PDE proposals and eight budget-matched space-filling proposals are frozen, with
one shared state and 15 unique proposed simulations. The 16 fixed-test outcomes
remain blind, and the saved probabilities are not yet claimed to be empirically
calibrated. See `docs/ML_BC_PROBABILISTIC_BOUNDARY_RESULTS.md`.

The completed budget-matched comparison shows why marginal classifier
uncertainty is not the final scientific objective. Both frozen GP models classify
all 16 blind test states correctly, while the uncertainty-selected model has the
lower Brier score (`0.01594` versus `0.02936`), lower log loss (`0.08682` versus
`0.12836`), and about nine times the reduction in the predefined
separator-band uncertainty per new PDE run. Nevertheless, none of its four
label-blind endpoint proposals is a B--C bracket; the space-filling control
finds one. That single control bracket dynamically recovers the independently
known `E_BC` orbit. The supported conclusion is therefore that reducing
marginal boundary uncertainty is not sufficient for constructing useful edge
brackets. See `docs/ML_BC_CONTROLLED_EXPERIMENT_RESULTS.md`.

`ml_BC_pair_acquisition_*_v1` replaces point acquisition with a joint-posterior
pair criterion. It estimates the probability that two nearby latent GP values
have opposite signs, then applies distance, local-normal orientation, and batch
diversity factors. Hyperparameter sensitivity across the saved model, a
wide-bound MLE, and a weakly regularized wide-bound MAP selects the same frozen
four-pair batch. With eight new endpoint integrations, three of four pairs are
true B--C brackets. All three complete ten binary bisections and recover the
independently Newton-converged edge orbit: minimum orbit distances are between
`6.66e-5` and `1.57e-4`, shadow durations are `146.5`--`162.5`, and direct
period errors are below `0.185%`. This is a retrospective B--C method-development
benchmark, not a held-out confirmation. The narrow `lambda approximately
0.276` C--B--A interval remains untouched for the frozen discovery test. See
`docs/ML_BC_PAIR_ACQUISITION_RESULTS.md`.

The frozen `ml_heldout_CBA_adaptive_v1` challenge then tests whether the method
can recover the narrow third-basin intrusion near `lambda = 0.276`. The
prespecified eight-run endpoint is negative: seven new states reach A, one
reaches C, and none reaches B. The final ML pair is nevertheless a C/A interval
of width `9.765625e-4` that contains the historically known B pilot. A separate
post-budget, five-way guarded refinement requires four more PDE runs to find B
at `lambda = 0.2761840820` and produces C--B and B--A brackets, each of width
`6.1035e-5`. The failure and follow-up are reported separately; the extra runs
do not count toward the frozen endpoint. The next dynamical tests are recovery
of known `E_BC` from the C--B bracket and unconstrained edge tracking of the
B--A bracket. See `docs/ML_HELDOUT_CBA_DISCOVERY_RESULTS.md`.

The post-budget C--B bracket now passes that dynamical test. A one-cycle pilot
initially remains far from `E_BC`, but three further full-state rebracketing
cycles produce edge-shadow durations `89.8`, `88.2`, and `92.0` and a minimum
physical distance `3.55e-5` to the independently Newton-converged periodic
edge orbit. Thus the narrow held-out-region C--B boundary is connected to the
stable manifold of the known `E_BC`. The nominal B--A continuation does not
remain binary: after one divergence, its evolved-state chord yields B, B, then
C and stops under the five-way guard. No direct B--A edge state is claimed;
that chord requires a new multiclass scan before further pairwise tracking. See
`docs/ML_HELDOUT_CBA_EDGE_TRACKING_RESULTS.md`.

`ml_heldout_CBA_evolved_BA_scan_v1` now performs that scan on the exact evolved
endpoint chord at global edge time `7.6`. Its 19 nonuniform samples give 12 B,
2 C, and 5 A outcomes, with no R5, R6, or final unresolved cases. The sampled
ordering is B through `alpha=0.825`, C at `0.85` and `0.875`, and A from `0.9`
through `1`. Thus the nominal evolved B--A chord contains a finite sampled C
region and has candidate B--C and C--A intervals of width `0.025`; these are
resolution-qualified candidates rather than proven direct adjacencies. The
recommended bounded closure is five-way refinement and a short dynamic pilot
of the C--A side, after which the project can pause without beginning another
exact-orbit continuation campaign. See `docs/EVOLVED_BA_MULTICLASS_SCAN_RESULTS.md`
and `docs/PROJECT_STATUS_AND_TEMPORARY_CONCLUSION.md`.

That bounded closure now terminates at its first guarded midpoint:
`alpha=0.8875` reaches B by frozen time `200`, with late median B distance
`2.5544e-3`. The apparent C--A interval therefore contains another B intrusion,
giving the sampled local ordering B--C--B--A. In accordance with the frozen
stopping rule, no dynamic C--A edge tracking is run. This result supports
multiclass folded/tongue-like basin intersections but does not establish
fractal, riddled, or Wada geometry. The present computational phase is now
closed; see `docs/BOUNDED_CA_CLOSURE_RESULTS.md`.

The small four-pair ML comparison is now followed by one separately versioned
robustness benchmark. `ml_BC_matched_benchmark_protocol_v1` freezes 12 matched,
outcome-blind randomized-maximin initial designs before completing the remaining
labels on the existing 17-by-17 B--C development plane. Each replicate gives
space filling, marginal pointwise uncertainty, and joint-posterior pair
acquisition the same 25-state training design and the same budget of four
disjoint pairs/eight newly revealed endpoint labels. Primary inference is at
the initial-design level; conditional guarded refinement and exact `E_BC`
recovery use only the highest-ranked valid pair per method/design. The complete
prespecified design is documented in `docs/ML_BC_MATCHED_BENCHMARK_DESIGN.md`.

The matched benchmark is complete. Across the 12 designs, pair-aware
acquisition produces 32 valid B--C brackets among 48 proposals, compared with
16 for pointwise uncertainty and four for space filling under identical
endpoint-label budgets. Pair-aware acquisition beats or ties both comparators
in every matched replicate. The 26 conditionally selected method/replicate
brackets reduce to 20 unique physical pairs; all 20 pass ten-step B/C-only
guarded refinement and recover the independently Newton-converged periodic edge
orbit. This is a controlled masked-label replay result on one fixed
two-dimensional B--C plane, not a full-state or rare-basin discovery claim.
See `docs/ML_BC_MATCHED_BENCHMARK_RESULTS.md` for the frozen interpretation,
statistics, integrity hashes, and manuscript claims.

The probability-index audit and corrected robustness extension supersede the
original `32/48`, `16/48`, `4/48` comparison for manuscript inference. Across
40 matched designs on the same masked development plane, corrected pair-aware
acquisition produces 112 valid brackets among 160 proposals, compared with 56
for corrected pointwise uncertainty and none for space filling. The paired
advantage persists across the tested initial-label budgets, kernel families,
and pair-distance scales. Ablation identifies joint opposite-sign probability
as the main task-aligned term, while the distance factor controls a
yield--length tradeoff and the orientation factor adds no detectable yield.
See `docs/ML_BC_MATCHED_BENCHMARK_PROBABILITY_INDEX_AUDIT.md` and
`docs/ML_BC_ACQUISITION_ROBUSTNESS_RESULTS.md`.

A final masked-replay consistency extension makes `qD`--the parsimonious score
frozen before prospective transfer--the primary rule across the same 40
designs and sensitivity axes. The `qD` reference produces 111 valid brackets
among 160 proposals, compared with 56 for pointwise uncertainty and none for
space filling; its paired advantage remains positive across every tested
initial-label count, kernel, and distance scale. The one-bracket difference
from `qDO` is not detectable at the replicate level, while unpenalized `q`
achieves greater raw yield by selecting substantially longer pairs. Manuscript
inference should therefore use `111/160` for the final `qD` method and retain
`q` and `qDO` as yield--length and orientation ablations. See
`docs/ML_BC_QD_PRIMARY_ROBUSTNESS_RESULTS.md`.

A separately frozen prospective transfer test now completes the end-to-end
workflow on a new physical-chord/POD3 state family rather than replaying the
original PC1--PC2 label bank. With the same eight-endpoint budget per method,
pair-aware `qD` produces two valid B--C brackets among four proposals;
corrected pointwise uncertainty and space filling produce none. Ten-step
five-way guarded refinement of the prespecified pair-aware bracket remains
B/C-only and recovers the independently computed `E_BC` neighborhood with
minimum phase-aware distance `1.3143e-4`, shadow duration `182.5`, and relative
period error `1.4221e-3`. This is one descriptive prospective transfer result
within a prelocalized pairwise family, not general or rare-basin discovery.
See `docs/ML_BC_PROSPECTIVE_TRANSFER_RESULTS.md`.
