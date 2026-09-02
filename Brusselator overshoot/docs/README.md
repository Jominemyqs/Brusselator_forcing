# Experiment documentation

This directory contains the versioned design decisions, validation records,
result summaries, and audit notes associated with the numerical experiments.
The repository-root `README.md` gives the chronological scientific overview and
links to the relevant documents here.

The manuscript and its literature synthesis are maintained separately in the
sibling `../../paper/` workspace so that this repository remains focused on
simulation, analysis, and reproducible experiment records.

## Repository-level guides

- `REPOSITORY_GUIDE.md`: stable source-tree map and principal workflows.
- `REPRODUCIBILITY.md`: software environment, data separation, and archive
  restoration policy.
- `MANUSCRIPT_DATA_MANIFEST.csv`: checksums for compact tables consumed by the
  manuscript Supplement.
- `CNSNS_RELEASE_CHECKLIST.md`: remaining code/data tasks before submission.

The other files in this directory are versioned experiment designs, results,
and implementation audits. Their filenames retain the experiment identifiers
used by the corresponding `experiment_outputs/` directories.
