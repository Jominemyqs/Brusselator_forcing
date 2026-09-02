#!/usr/bin/env python3
"""Validate the frozen, budget-matched B--C comparison artifacts."""

import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parent
OUTDIR = ROOT / "experiment_outputs" / "ml_BC_budget_comparison_v1"


def sha256_file(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main():
    manifest = json.loads((OUTDIR / "frozen_model_manifest.json").read_text(encoding="utf-8"))
    summary = json.loads((OUTDIR / "model_summary.json").read_text(encoding="utf-8"))
    config = json.loads((OUTDIR / "ml_bc_budget_comparison_config.json").read_text(encoding="utf-8"))
    active_training = pd.read_csv(OUTDIR / "active_training_set.csv")
    baseline_training = pd.read_csv(OUTDIR / "baseline_training_set.csv")
    active_blind = pd.read_csv(OUTDIR / "active_fixed_test_blind_predictions.csv")
    baseline_blind = pd.read_csv(OUTDIR / "baseline_fixed_test_blind_predictions.csv")
    source_geometry = pd.read_csv(ROOT / config["geometry_manifest_csv"])
    active_proposals = pd.read_csv(ROOT / config["active_proposals_csv"])
    baseline_proposals = pd.read_csv(ROOT / config["baseline_proposals_csv"])
    active_model = np.load(OUTDIR / "active_probit_gp_model.npz")
    baseline_model = np.load(OUTDIR / "baseline_probit_gp_model.npz")
    active_brackets = pd.read_csv(OUTDIR / "active_edge_bracket_proposals.csv")
    baseline_brackets = pd.read_csv(OUTDIR / "baseline_edge_bracket_proposals.csv")
    label_like = {"final_label", "outcome_label", "true_label", "predicted_label"}
    active_seed = set(active_training.loc[active_training.design_role == "initial_seed", "state_index"])
    baseline_seed = set(baseline_training.loc[baseline_training.design_role == "initial_seed", "state_index"])
    active_added = set(active_training.state_index) - active_seed
    baseline_added = set(baseline_training.state_index) - baseline_seed
    expected_test = list(source_geometry.loc[source_geometry.design_role == "fixed_test", "state_index"])

    checks = {
        "models_declared_frozen": manifest["models_frozen"] is True,
        "fixed_test_labels_not_used": manifest["fixed_test_labels_used"] is False,
        "config_has_no_test_label_input": not any("test_label" in key for key in config),
        "active_count_is_33": len(active_training) == 33 == active_model["X"].shape[0],
        "baseline_count_is_33": len(baseline_training) == 33 == baseline_model["X"].shape[0],
        "shared_seed_is_same_25": len(active_seed) == 25 and active_seed == baseline_seed,
        "active_additions_are_frozen_eight": active_added == set(active_proposals.state_index),
        "baseline_additions_are_frozen_eight": baseline_added == set(baseline_proposals.state_index),
        "proposal_overlap_is_one": len(active_added & baseline_added) == 1,
        "training_labels_are_binary": set(active_training.final_label) <= {"B", "C"}
        and set(baseline_training.final_label) <= {"B", "C"},
        "blind_prediction_counts_are_16": len(active_blind) == len(baseline_blind) == 16,
        "blind_state_order_matches": list(active_blind.state_index) == list(baseline_blind.state_index) == expected_test,
        "blind_files_contain_no_labels": not (label_like & set(active_blind.columns))
        and not (label_like & set(baseline_blind.columns)),
        "blind_probabilities_are_valid": active_blind.probability_C.between(0, 1).all()
        and baseline_blind.probability_C.between(0, 1).all(),
        "active_model_hash_matches": sha256_file(OUTDIR / "active_probit_gp_model.npz")
        == manifest["active_model_hash"],
        "baseline_model_hash_matches": sha256_file(OUTDIR / "baseline_probit_gp_model.npz")
        == manifest["baseline_model_hash"],
        "active_prediction_hash_matches": sha256_file(OUTDIR / "active_fixed_test_blind_predictions.csv")
        == manifest["active_blind_prediction_hash"],
        "baseline_prediction_hash_matches": sha256_file(OUTDIR / "baseline_fixed_test_blind_predictions.csv")
        == manifest["baseline_blind_prediction_hash"],
        "active_bracket_hash_matches": sha256_file(OUTDIR / "active_edge_bracket_proposals.csv")
        == manifest["active_edge_bracket_proposal_hash"],
        "baseline_bracket_hash_matches": sha256_file(OUTDIR / "baseline_edge_bracket_proposals.csv")
        == manifest["baseline_edge_bracket_proposal_hash"],
        "four_label_blind_brackets_per_model": len(active_brackets) == len(baseline_brackets) == 4,
        "active_brackets_straddle_half": (active_brackets.B_probability_C < 0.5).all()
        and (active_brackets.C_probability_C >= 0.5).all(),
        "baseline_brackets_straddle_half": (baseline_brackets.B_probability_C < 0.5).all()
        and (baseline_brackets.C_probability_C >= 0.5).all(),
        "summary_counts_match": summary["active_training_count"] == summary["baseline_training_count"] == 33,
        "exact_edge_excluded_from_fit": summary["exact_edge_orbit_used_for_fitting"] is False,
        "heldout_CA_excluded": summary["held_out_CA_region_used"] is False,
    }
    failed = [name for name, passed in checks.items() if not passed]
    pd.DataFrame([{"check": name, "passed": passed} for name, passed in checks.items()]).to_csv(
        OUTDIR / "validation_checks.csv", index=False
    )
    if failed:
        raise RuntimeError("Budget-comparison validation failed: " + ", ".join(failed))
    print(f"All {len(checks)} budget-comparison checks passed.")


if __name__ == "__main__":
    main()
