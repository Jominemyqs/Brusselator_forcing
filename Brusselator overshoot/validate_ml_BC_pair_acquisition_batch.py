#!/usr/bin/env python3
"""Validate the outcome-blind, frozen bracket-aware development batch."""

import hashlib
import json
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parent
OUTDIR = ROOT / "experiment_outputs" / "ml_BC_pair_acquisition_frozen_v1"


def sha256_file(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main():
    manifest = json.loads((OUTDIR / "frozen_pair_manifest.json").read_text(encoding="utf-8"))
    pairs = pd.read_csv(OUTDIR / "frozen_pair_batch.csv")
    endpoints = pd.read_csv(OUTDIR / "frozen_endpoint_manifest.csv")
    geometry = pd.read_csv(
        ROOT / "experiment_outputs" / "ml_BC_boundary_geometry_v1" / "sampling_state_manifest.csv"
    )
    role = geometry.set_index("state_index").loc[endpoints.state_index, "design_role"]
    forbidden = {"final_label", "outcome_label", "true_label", "raw_five_way_outcome"}
    checks = {
        "batch_declared_frozen": manifest["pair_batch_frozen"] is True,
        "endpoint_outcomes_unused": manifest["endpoint_outcomes_used"] is False,
        "four_pairs": len(pairs) == 4 and pairs.pair_id.nunique() == 4,
        "eight_unique_endpoints": len(endpoints) == endpoints.state_index.nunique() == 8,
        "endpoint_budget_is_eight": manifest["new_PDE_endpoint_budget"] == 8,
        "all_endpoints_are_acquisition_pool": (role == "acquisition_pool").all(),
        "batch_contains_no_labels": not (forbidden & set(pairs.columns)),
        "endpoint_manifest_contains_no_labels": not (forbidden & set(endpoints.columns)),
        "pair_hash_matches": sha256_file(OUTDIR / "frozen_pair_batch.csv")
        == manifest["frozen_pair_batch_hash"],
        "endpoint_hash_matches": sha256_file(OUTDIR / "frozen_endpoint_manifest.csv")
        == manifest["frozen_endpoint_manifest_hash"],
        "ranking_robustness_gate_passed": manifest["minimum_top_20_jaccard_across_sensitivity_fits"] >= 0.8
        and manifest["four_pair_batch_intersection_across_every_fit"] == 4,
        "exact_edge_excluded": manifest["exact_edge_orbit_used"] is False,
        "heldout_CA_excluded": manifest["held_out_CA_region_used"] is False,
    }
    failed = [name for name, passed in checks.items() if not passed]
    if failed:
        raise RuntimeError("Frozen pair-batch validation failed: " + ", ".join(failed))
    print(f"All {len(checks)} frozen pair-batch checks passed.")


if __name__ == "__main__":
    main()
