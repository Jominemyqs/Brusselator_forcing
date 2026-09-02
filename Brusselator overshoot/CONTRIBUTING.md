# Contributing

This is a research-code repository whose versioned configurations and result
notes support a scientific manuscript. Changes should preserve numerical
provenance and keep scientific definitions explicit.

## Before changing code

1. Read `README.md`, `docs/REPOSITORY_GUIDE.md`, and the relevant experiment
   design/result note.
2. Check `git status`; do not overwrite an existing experiment directory.
3. Treat the saved configuration and machine-readable summaries as the source
   of truth for a completed experiment.
4. Do not silently change physical norms, phase/reflection handling,
   classification thresholds, random seeds, solver tolerances, or outcome
   definitions.

## New experiments

- Create a new versioned directory below `experiment_outputs/`.
- Save configuration, runtime metadata, dependency/provenance information,
  compact summaries, and raw output separately.
- Make the runner refuse to overwrite an existing output directory.
- Record inconclusive or failed cases rather than removing them.
- Add a concise design or results note to `docs/` when the experiment affects
  a manuscript claim.

## Validation

Run:

```bash
python -m pytest
python scripts/verify_repository.py
```

For MATLAB changes, also run the smallest relevant configuration or numerical
cross-check. Expensive experiments should not be repeated merely for code-style
changes.

## Generated data

Raw `.mat` files and versioned output directories are ignored by default. Do
not force-add large files to Git. Select publication data should instead be
archived with a persistent identifier and checked against
`docs/MANUSCRIPT_DATA_MANIFEST.csv`.
