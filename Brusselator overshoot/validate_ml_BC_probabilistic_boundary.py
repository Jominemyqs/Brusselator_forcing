#!/usr/bin/env python3
"""Read-only validation of the first probabilistic B--C boundary output."""

import json
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parent


def main():
    with (ROOT / "ml_bc_probabilistic_config.json").open("r", encoding="utf-8") as handle:
        config = json.load(handle)
    outdir = ROOT / config["output_root"] / config["experiment_name"]
    required = [
        "model_summary.json",
        "probit_gp_model.npz",
        "training_posterior_predictions.csv",
        "nested_LOO_predictions.csv",
        "acquisition_pool_predictions.csv",
        "active_acquisition_proposals.csv",
        "space_filling_baseline_proposals.csv",
        "fixed_test_blind_predictions.csv",
        "dense_probability_grid.csv",
        "probability_half_contour.csv",
        "probabilistic_boundary_and_acquisition.png",
    ]
    with (outdir / "model_summary.json").open("r", encoding="utf-8") as handle:
        summary = json.load(handle)
    with (outdir / "ml_bc_probabilistic_config.json").open("r", encoding="utf-8") as handle:
        saved_config = json.load(handle)
    training = pd.read_csv(outdir / "training_posterior_predictions.csv")
    loo = pd.read_csv(outdir / "nested_LOO_predictions.csv")
    pool = pd.read_csv(outdir / "acquisition_pool_predictions.csv")
    active = pd.read_csv(outdir / "active_acquisition_proposals.csv")
    baseline = pd.read_csv(outdir / "space_filling_baseline_proposals.csv")
    test = pd.read_csv(outdir / "fixed_test_blind_predictions.csv")
    contour = pd.read_csv(outdir / "probability_half_contour.csv")
    model = np.load(outdir / "probit_gp_model.npz")

    checks = {
        "all_required_outputs_exist": all((outdir / name).is_file() for name in required),
        "training_count_is_25": len(training) == 25 and len(loo) == 25,
        "training_is_only_B_C": set(training["final_label"]) == {"B", "C"},
        "pool_count_is_248": len(pool) == 248,
        "fixed_test_count_is_16": len(test) == 16,
        "saved_configuration_matches": saved_config == config,
        "fixed_test_has_no_labels": not any(
            name in test.columns for name in ["final_label", "raw_five_way_outcome", "outcome_label"]
        ),
        "fixed_test_is_still_unrun": np.all(test["simulation_status"] == "not_run"),
        "active_batch_has_requested_size": len(active) == config["active_batch_size"],
        "baseline_has_requested_size": len(baseline) == config["active_batch_size"],
        "proposal_indices_are_unique": active["state_index"].is_unique
        and baseline["state_index"].is_unique,
        "proposals_come_only_from_pool": set(active["state_index"]).issubset(set(pool["state_index"]))
        and set(baseline["state_index"]).issubset(set(pool["state_index"])),
        "probabilities_are_valid": all(
            np.all((frame["probability_C"] >= 0.0) & (frame["probability_C"] <= 1.0))
            for frame in [training, pool, test]
        ),
        "entropies_are_valid": all(
            np.all((frame["predictive_entropy_bits"] >= 0.0) & (frame["predictive_entropy_bits"] <= 1.0 + 1e-12))
            for frame in [training, pool, test]
        ),
        "contour_is_nonempty": len(contour) > 1 and contour["segment_id"].nunique() >= 1,
        "model_dimensions_match": model["X"].shape == (25, 2)
        and model["y"].shape == (25,)
        and model["lengthscales"].shape == (2,),
        "laplace_converged": bool(summary["laplace_converged"]),
        "full_model_hyperparameters_are_interior": not bool(summary["signal_std_at_upper_bound"])
        and not bool(summary["lengthscale_at_bound"]),
        "fixed_test_labels_not_used": not bool(summary["fixed_test_labels_used"]),
        "exact_edge_excluded": not bool(summary["exact_edge_orbit_used"]),
        "held_out_region_excluded": not bool(summary["held_out_CA_region_used"]),
        "no_random_seed_needed": not bool(summary["random_seed_used"]),
    }
    failed = [name for name, value in checks.items() if not value]
    if failed:
        raise AssertionError("Probabilistic boundary checks failed: " + ", ".join(failed))
    for name, value in checks.items():
        print(f"{name}: {value}")
    print(f"All {len(checks)} probabilistic-boundary checks passed.")


if __name__ == "__main__":
    main()
