#!/usr/bin/env python3
"""Generate hash-frozen qD-primary masked-replay acquisition batches."""

from pathlib import Path

import run_ml_BC_acquisition_robustness_proposals as implementation


ROOT = Path(__file__).resolve().parent
implementation.CONFIG_FILE = ROOT / "ml_bc_acquisition_qd_primary_config.json"


if __name__ == "__main__":
    implementation.main()
