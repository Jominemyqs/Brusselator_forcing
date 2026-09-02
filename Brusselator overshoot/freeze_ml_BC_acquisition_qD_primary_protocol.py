#!/usr/bin/env python3
"""Freeze the qD-primary replay using the audited robustness machinery."""

from pathlib import Path

import freeze_ml_BC_acquisition_robustness_protocol as implementation


ROOT = Path(__file__).resolve().parent
implementation.CONFIG_FILE = ROOT / "ml_bc_acquisition_qd_primary_config.json"
implementation.IMPLEMENTATION_FILES = (
    "freeze_ml_BC_acquisition_qD_primary_protocol.py",
    "run_ml_BC_acquisition_qD_primary_proposals.py",
    "evaluate_ml_BC_acquisition_qD_primary.py",
    "freeze_ml_BC_acquisition_robustness_protocol.py",
    "run_ml_BC_acquisition_robustness_proposals.py",
    "evaluate_ml_BC_acquisition_robustness.py",
    "probit_gp_kernels.py",
)


if __name__ == "__main__":
    implementation.main()
