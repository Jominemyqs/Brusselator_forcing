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

The Git repository contains code, configuration files, and human-readable
experiment records. The local `experiment_outputs/` tree contains several
gigabytes of generated trajectories and is intentionally ignored by Git.

`docs/MANUSCRIPT_DATA_MANIFEST.csv` lists the compact machine-readable result
files used to generate the manuscript Supplement, together with SHA-256 hashes
and byte counts. `scripts/verify_repository.py` checks the files when the local
data tree is available and reports missing files as an archive-restoration
requirement rather than a code failure.

## Public release gate

Before journal submission or public release:

1. deposit the selected raw and processed numerical data in Zenodo or another
   durable research-data archive;
2. create a versioned Git release for the exact code snapshot;
3. record both persistent identifiers in `CITATION.cff`, this file, the
   manuscript data-availability statement, and the Supplement;
4. verify the restored archive against `docs/MANUSCRIPT_DATA_MANIFEST.csv`;
5. choose and add an explicit software license.

No DOI or license is asserted in the working repository until those choices
have actually been completed.

## Scope of reproducibility

The archived calculations support fixed-domain semidiscrete numerical claims
for the declared grids, solvers, tolerances, seeds, and finite horizons. They do
not constitute a continuum existence theorem, a proof of global basin
structure, or a claim that the tested forcing family exhausts all physically
possible protocols.
