# Generated experiment outputs

Each subdirectory is a versioned numerical study produced by the scripts in the
repository root. The files remain the local source of truth for the reported
results, but the directory is ignored by Git because it contains several
gigabytes of generated MATLAB states and derived artifacts.

Each new experiment should continue to save its configuration, runtime
metadata, compact summaries, and raw numerical data separately. A deliberately
selected small artifact can be added with `git add -f`, but large raw output
should normally be archived outside ordinary Git history.

The compact files directly consumed by the manuscript Supplement are listed in
`../docs/MANUSCRIPT_DATA_MANIFEST.csv`. Before public release, selected raw and
processed outputs should be deposited in a persistent research-data archive;
the archive must be tested against that manifest before its DOI is cited.
