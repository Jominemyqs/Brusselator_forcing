#!/usr/bin/env python3
"""Freeze the physical-chord transfer protocol before any new PDE labels."""

from __future__ import annotations

import hashlib
import json
import shutil
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parent
CONFIG = ROOT / "ml_bc_transfer_prospective_config.json"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    config = json.loads(CONFIG.read_text(encoding="utf-8"))
    outdir = ROOT / config["protocol_directory"]
    if outdir.exists():
        raise FileExistsError(f"Refusing to overwrite frozen transfer protocol: {outdir}")
    geometry_csv = ROOT / config["geometry_manifest_csv"]
    geometry_mat = ROOT / config["geometry_mat"]
    geometry_summary = ROOT / config["geometry_directory"] / "transfer_geometry_summary.csv"
    geometry = pd.read_csv(geometry_csv)
    summary = pd.read_csv(geometry_summary).iloc[0]
    gate = config["geometry_gate"]
    if len(geometry) != gate["expected_state_count"]:
        raise ValueError("Transfer state count changed")
    if int((geometry.design_role == "initial_seed").sum()) != gate["expected_initial_count"]:
        raise ValueError("Transfer initial design count changed")
    if geometry.outcome_label.fillna("").astype(str).str.len().sum() != 0:
        raise ValueError("A transfer candidate label was populated before protocol freeze")
    if set(geometry.simulation_status) != {"not_run"}:
        raise ValueError("A transfer candidate was simulated before protocol freeze")
    if int(summary.B_anchor_state_index) != gate["expected_B_anchor_state_index"] or int(
        summary.C_anchor_state_index
    ) != gate["expected_C_anchor_state_index"]:
        raise ValueError("Physical-chord anchor indices changed")
    if float(summary.projector_difference_from_PC1_PC2) < gate["minimum_projector_difference_from_PC1_PC2"]:
        raise ValueError("Transfer family is insufficiently distinct from PC1--PC2")

    initial = geometry[geometry.design_role == "initial_seed"].copy()
    initial.insert(0, "selection_position", range(1, len(initial) + 1))
    outdir.mkdir(parents=True)
    shutil.copy2(CONFIG, outdir / CONFIG.name)
    initial_path = outdir / "frozen_initial_state_manifest.csv"
    initial.to_csv(initial_path, index=False)
    manifest = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "protocol_frozen": True,
        "protocol_version": "ml_BC_transfer_protocol_v1",
        "geometry_id": gate["expected_geometry_id"],
        "geometry_manifest_sha256": sha256(geometry_csv),
        "geometry_mat_sha256": sha256(geometry_mat),
        "geometry_summary_sha256": sha256(geometry_summary),
        "configuration_sha256": sha256(CONFIG),
        "initial_state_manifest_sha256": sha256(initial_path),
        "initial_state_count": len(initial),
        "pair_count_per_method": config["batch"]["pair_count"],
        "endpoint_budget_per_method": config["batch"]["endpoint_budget_per_method"],
        "candidate_outcomes_used": False,
        "nonanchor_outcomes_used_for_geometry": False,
        "source_anchor_labels_used_only_to_prelocalize_family": True,
        "exact_edge_orbit_used": False,
        "development_label_bank_used": False,
        "heldout_CBA_region_used": False,
        "implementation_sha256": {Path(__file__).name: sha256(Path(__file__))},
        "interpretation": "frozen prospective transfer on one physical-chord/POD3 state family",
    }
    manifest_path = outdir / "frozen_transfer_protocol_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(json.dumps(manifest, indent=2))
    print(f"Frozen prospective transfer protocol saved in: {outdir}")


if __name__ == "__main__":
    main()
