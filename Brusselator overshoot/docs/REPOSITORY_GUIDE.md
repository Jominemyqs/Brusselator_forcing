# Repository guide

## Why the source files remain at repository root

The code is organized logically rather than by moving the existing MATLAB and
Python files into new source directories. Completed experiments use frozen
repository-root-relative paths such as `experiment_outputs/...`; moving the
entry points would change those assumptions and make the recorded provenance
harder to audit. New organizational material therefore indexes the stable
source tree without rewriting it.

## Source categories

| Pattern | Role |
| --- | --- |
| `brusselator_*.m` | Shared model, configuration, norms, classification, POD, and edge-tracking utilities |
| `solve_brusselator_1d_forced.m` | Forced PDE integration |
| `run_*.m` | MATLAB experiment entry points |
| `analyze_*.m`, `build_*.m` | MATLAB postprocessing and geometry construction |
| `freeze_*.py` | Freeze acquisition protocols or proposal batches before labels are opened |
| `run_*.py`, `fit_*.py` | Python acquisition/model execution |
| `evaluate_*.py`, `analyze_*.py`, `validate_*.py` | Python evaluation, audits, and validation |
| `probit_gp*.py`, `multiclass_probit_gp.py` | Probabilistic classifier implementations |
| `*_config.m`, `*_config.json` | Versioned scientific and computational definitions |
| `docs/` | Human-readable designs, results, audits, and release documentation |
| `experiment_outputs/` | Local generated data; excluded from ordinary Git history |

Legacy course-era scripts and top-level output directories remain in place for
historical reproducibility. They are not the preferred entry points for the
current manuscript.

## Principal manuscript workflow

### 1. Frozen outcomes and forcing landings

1. `brusselator_default_config.m`
2. `brusselator_execute_forcing_protocol.m`
3. `run_candidate_c_numerical_validation.m`
4. `run_recurrent_class_validation.m`
5. `build_stage2_landing_geometry.m`

### 2. Multiclass slices and the periodic B/C edge state

1. `run_five_way_edge_slice_scan.m`
2. `run_edge_slice_pairwise_refinement.m`
3. `run_dynamic_edge_tracking_BC_multicycle.m`
4. `run_periodic_edge_orbit_newton.m`
5. `run_periodic_edge_orbit_floquet.m`
6. `run_periodic_edge_orbit_unstable_branches.m`
7. `run_BC_edge_resolution_geometry.m`
8. `run_physical_ramp_chord_BC_edge.m`

### 3. Development-family acquisition benchmark

1. `build_ml_BC_boundary_geometry.m`
2. `run_ml_BC_boundary_initial_labeling.m`
3. `freeze_ml_BC_acquisition_qD_primary_protocol.py`
4. `run_ml_BC_acquisition_qD_primary_proposals.py`
5. `evaluate_ml_BC_acquisition_qD_primary.py`

The complete 17-by-17 label bank is a local generated artifact and must be
restored from the publication data archive before replaying all 40 matched
designs.

### 4. Prospective transfer and held-out challenge

1. `build_ml_BC_transfer_geometry.m`
2. `freeze_ml_BC_transfer_protocol.py`
3. `run_ml_BC_transfer_initial_labeling.m`
4. `freeze_ml_BC_transfer_proposals.py`
5. `run_ml_BC_transfer_endpoint_labeling.m`
6. `evaluate_ml_BC_transfer_prospective.py`
7. `run_ml_BC_transfer_edge_recovery.m`

The separate rare-basin challenge is defined by
`freeze_ml_heldout_CBA_protocol.py` and the corresponding
`run_ml_heldout_CBA_*` entry points. Its frozen-budget failure is retained as a
scope limitation, not tuned away.

## Scientific terminology

- An **opposite-label pair** has endpoints classified as different outcomes.
- A **guarded bracket** is an opposite-label pair whose refinement remains
  pairwise under the full five-outcome classifier.
- `E_BC` denotes the Newton-converged periodic edge orbit of the stated
  semidiscrete problem.
- POD coordinates are used for visualization and state-family construction;
  physical recurrence and classification use the declared quadrature-weighted
  full-state norm.
