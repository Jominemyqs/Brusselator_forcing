# CNSNS code-and-data release checklist

This checklist covers the research-software side of a submission to
*Communications in Nonlinear Science and Numerical Simulation*. Journal portal
fields and author declarations are tracked in the separate manuscript
workspace.

## Repository

- [x] Initialize Git on branch `main`.
- [x] Exclude raw numerical outputs, caches, editor backups, and local
  environments.
- [x] Add pinned Python dependencies and a lightweight verification command.
- [x] Add code-navigation and reproducibility documentation.
- [x] Add citation metadata without inventing a DOI.
- [ ] Select a software license with the author/advisor.
- [x] Add the configured GitHub repository URL to `CITATION.cff`.
- [ ] Tag the submitted snapshot (suggested tag: `cnsns-submission-v1`).

## Data archive

- [ ] Decide which raw trajectories and processed tables are required for
  reproduction versus audit only.
- [x] Bundle all 18 compact supplementary-table sources outside the raw-output tree.
- [ ] Bundle additional raw trajectories and frozen inputs for full reproduction.
- [x] Verify all compact manuscript sources against
  `MANUSCRIPT_DATA_MANIFEST.csv`.
- [ ] Deposit the archive and record its DOI.
- [ ] Test restoration in a clean directory.

## Manuscript linkage

- [x] Add the repository citation and exact code revision to the manuscript and Supplement.
- [ ] Add an archival data DOI when deposited.
- [x] Confirm that repository terminology matches the paper: opposite-label
  pair, guarded bracket, signed unstable-direction perturbation, and `E_BC`.
- [x] Confirm that the paper reports `qD` as the final method and `qDO` only as
  an ablation.
- [ ] Freeze the exact manuscript, Supplement, figure, code, and data hashes
  used at submission.
