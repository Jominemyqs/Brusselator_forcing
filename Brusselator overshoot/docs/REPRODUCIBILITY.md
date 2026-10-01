# Reproducibility and data availability

## Software environment

The numerical study uses MATLAB for the reaction--diffusion integrations,
shooting, Floquet calculations, and state-space geometry, and Python for the
probabilistic classifier and acquisition benchmarks.

- MATLAB: the final repository audit was performed with R2026a. Individual
  experiment directories retain the runtime metadata of the original run.
- Python: the release audit passes under Python 3.9.6. Exact analysis
  dependencies are pinned in `requirements.txt`; test dependencies are in
  `requirements-dev.txt`.

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -r requirements-dev.txt
python -m pytest
python scripts/verify_repository.py
```

MATLAB entry points must be run from the repository root because completed
experiments use frozen paths relative to that directory.

## Code, compact evidence, and raw data

Repository: https://github.com/Jominemyqs/Brusselator_forcing

Version 0.1.1 includes the post hoc pointwise-orientation ablation script and
`data/manuscript_compact_results.zip`. The ZIP contains the 18 compact result
files used to generate the current supplementary tables, byte-for-byte copies
of the frozen experiment outputs, and their checksum manifest. Restore it from
the repository root with:

```bash
python -m zipfile -e data/manuscript_compact_results.zip .
python scripts/verify_repository.py
```

`docs/MANUSCRIPT_DATA_MANIFEST.csv` records SHA-256 hashes and byte counts.
`data/README.txt` also records the ZIP checksum. The table generator is in the
separate manuscript workspace at `paper/supplement/scripts/build_tables.py`.
The restored compact files support table regeneration; they are not sufficient
alone to rerun all simulations, regenerate all figures, or replay every
acquisition experiment. Those operations additionally require the frozen
geometry, label banks, protocol records, and/or raw trajectories described by
the experiment configurations.

The several-gigabyte raw-output tree remains outside Git. No research-data
archive DOI or software license has been assigned. These limitations must not
be represented as a complete raw-data release.

## Before submission

1. Synchronize the versioned code and compact-data bundle with the public
   repository, and cite the exact Git revision in Supplement S1.3.
2. Select and deposit the additional data needed for full experimental
   reproduction in a durable archive; record its persistent identifier in the
   manuscript and repository when available. If data cannot be shared, state
   the reason in accordance with the journal's data policy.
3. Choose an explicit software license before claiming licensed reuse.
4. Tag the exact submitted snapshot after the planned scientific revisions.

## Scope of reproducibility

The archived calculations support fixed-domain semidiscrete numerical claims
for the declared grids, solvers, tolerances, seeds, and finite horizons. They do
not constitute a continuum existence theorem, a proof of global basin
structure, or a claim that the tested forcing family exhausts all physically
possible protocols.
