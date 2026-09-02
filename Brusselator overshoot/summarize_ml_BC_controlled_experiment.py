#!/usr/bin/env python3
"""Consolidate the four controlled B--C ML evaluation layers."""

import hashlib
import json
import os
from datetime import datetime, timezone
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/private/tmp/brusselator_matplotlib")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parent
MODEL_DIR = ROOT / "experiment_outputs" / "ml_BC_budget_comparison_v1"
EVALUATION_DIR = ROOT / "experiment_outputs" / "ml_BC_budget_evaluation_v1"
ENDPOINT_DIR = ROOT / "experiment_outputs" / "ml_BC_edge_bracket_labeling_v1"
EDGE_DIR = ROOT / "experiment_outputs" / "ml_BC_edge_recovery_v1"
OUTDIR = ROOT / "experiment_outputs" / "ml_BC_controlled_experiment_v1"


def sha256_file(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def frozen_gate():
    frozen = json.loads((MODEL_DIR / "frozen_model_manifest.json").read_text(encoding="utf-8"))
    checks = {
        "active_model": (MODEL_DIR / "active_probit_gp_model.npz", frozen["active_model_hash"]),
        "baseline_model": (MODEL_DIR / "baseline_probit_gp_model.npz", frozen["baseline_model_hash"]),
        "active_predictions": (
            MODEL_DIR / "active_fixed_test_blind_predictions.csv",
            frozen["active_blind_prediction_hash"],
        ),
        "baseline_predictions": (
            MODEL_DIR / "baseline_fixed_test_blind_predictions.csv",
            frozen["baseline_blind_prediction_hash"],
        ),
        "active_brackets": (
            MODEL_DIR / "active_edge_bracket_proposals.csv",
            frozen["active_edge_bracket_proposal_hash"],
        ),
        "baseline_brackets": (
            MODEL_DIR / "baseline_edge_bracket_proposals.csv",
            frozen["baseline_edge_bracket_proposal_hash"],
        ),
    }
    failures = [name for name, (path, expected) in checks.items() if sha256_file(path) != expected]
    if failures:
        raise RuntimeError("Frozen artifact changed: " + ", ".join(failures))
    return frozen


def bracket_table(endpoint_labels):
    label_map = endpoint_labels.set_index("state_index")["final_label"].to_dict()
    rows = []
    for model in ["active", "baseline"]:
        proposals = pd.read_csv(MODEL_DIR / f"{model}_edge_bracket_proposals.csv")
        for proposal in proposals.itertuples(index=False):
            first = label_map[int(proposal.B_state_index)]
            second = label_map[int(proposal.C_state_index)]
            verified = {first, second} == {"B", "C"}
            rows.append(
                {
                    "model": model,
                    "proposal_rank": int(proposal.proposal_rank),
                    "predicted_B_state_index": int(proposal.B_state_index),
                    "predicted_C_state_index": int(proposal.C_state_index),
                    "predicted_B_actual_label": first,
                    "predicted_C_actual_label": second,
                    "verified_BC_bracket": verified,
                    "predicted_orientation_correct": first == "B" and second == "C",
                }
            )
    return pd.DataFrame(rows)


def make_figure(path, cases, metrics, uncertainty, brackets):
    figure, axes = plt.subplots(1, 3, figsize=(15.5, 4.6), constrained_layout=True)
    x = np.arange(len(cases))
    target = (cases.true_label == "C").astype(float)
    axes[0].scatter(x, target, marker="s", s=48, c="black", label="PDE outcome")
    axes[0].plot(x, cases.active_probability_C, "o-", color="tab:green", label="active")
    axes[0].plot(x, cases.baseline_probability_C, "o-", color="tab:cyan", label="space-filling")
    axes[0].axhline(0.5, color="0.4", linestyle="--", linewidth=1)
    axes[0].set_xlabel("blind test case (state-index order)")
    axes[0].set_ylabel(r"$p_C$ / binary outcome")
    axes[0].set_title("Blind predictions: both 16/16")
    axes[0].legend(fontsize=8)

    model_names = ["active", "space-filling"]
    brier = [metrics["active"]["Brier_score_on_BC_cases"], metrics["space_filling_baseline"]["Brier_score_on_BC_cases"]]
    log_loss = [metrics["active"]["log_loss_on_BC_cases"], metrics["space_filling_baseline"]["log_loss_on_BC_cases"]]
    width = 0.36
    axes[1].bar(np.arange(2) - width / 2, brier, width, label="Brier")
    axes[1].bar(np.arange(2) + width / 2, log_loss, width, label="log loss")
    axes[1].set_xticks(np.arange(2), model_names)
    axes[1].set_ylabel("proper score (lower is better)")
    axes[1].set_title("Probability quality")
    axes[1].set_ylim(0, 0.16)
    axes[1].legend()
    for index, name in enumerate(["active", "baseline"]):
        reduction = uncertainty[name]["frozen_separator_band_normalized_integral_reduction_per_added_trajectory"]
        axes[1].text(index, max(brier[index], log_loss[index]) + 0.006, f"band Δ/run={reduction:.4g}", ha="center", fontsize=8)

    matrix = np.zeros((2, 4))
    for row, model in enumerate(["active", "baseline"]):
        selected = brackets[brackets.model == model].sort_values("proposal_rank")
        matrix[row, :] = selected.verified_BC_bracket.astype(float)
    axes[2].imshow(matrix, cmap="RdYlGn", vmin=0, vmax=1, aspect="auto")
    axes[2].set_xticks(np.arange(4), ["rank 1", "rank 2", "rank 3", "rank 4"])
    axes[2].set_yticks(np.arange(2), ["active", "space-filling"])
    axes[2].set_title("Independent PDE bracket check")
    for row in range(2):
        for column in range(4):
            axes[2].text(column, row, "B/C" if matrix[row, column] else "same", ha="center", va="center")
    axes[2].text(
        3, 1.37, "rank 4 → E_BC" + "\n" + "d_min = 1.74 × 10⁻⁴",
        ha="center", va="center", color="tab:blue", fontweight="bold"
    )
    figure.suptitle("Controlled, budget-matched B--C boundary experiment")
    figure.savefig(path, dpi=300)
    plt.close(figure)


def main():
    if OUTDIR.exists():
        raise FileExistsError(f"Refusing to overwrite consolidated output: {OUTDIR}")
    frozen = frozen_gate()
    model_summary = json.loads((MODEL_DIR / "model_summary.json").read_text(encoding="utf-8"))
    metrics = json.loads((EVALUATION_DIR / "blind_test_metrics.json").read_text(encoding="utf-8"))
    uncertainty = model_summary["uncertainty"]
    cases = pd.read_csv(EVALUATION_DIR / "blind_test_individual_predictions.csv")
    endpoint_labels = pd.read_csv(ENDPOINT_DIR / "edge_bracket_endpoint_labels.csv")
    edge = pd.read_csv(EDGE_DIR / "edge_recovery_summary.csv")
    brackets = bracket_table(endpoint_labels)

    rows = []
    for model, metric_key in [("active", "active"), ("baseline", "space_filling_baseline")]:
        edge_row = edge[edge.model == model].iloc[0]
        model_brackets = brackets[brackets.model == model]
        verified = model_brackets[model_brackets.verified_BC_bracket]
        hyper = model_summary[f"{model}_hyperparameters"]
        rows.append(
            {
                "model": model,
                "training_trajectory_budget": 33,
                "added_trajectory_budget": 8,
                "blind_correct_count": metrics[metric_key]["correct_count"],
                "blind_test_count": metrics[metric_key]["test_count"],
                "blind_accuracy": metrics[metric_key]["accuracy_all_16"],
                "blind_Brier_score": metrics[metric_key]["Brier_score_on_BC_cases"],
                "blind_log_loss": metrics[metric_key]["log_loss_on_BC_cases"],
                "full_plane_uncertainty_reduction_per_added_trajectory": uncertainty[model]["full_plane_normalized_integral_reduction_per_added_trajectory"],
                "separator_band_uncertainty_reduction_per_added_trajectory": uncertainty[model]["frozen_separator_band_normalized_integral_reduction_per_added_trajectory"],
                "verified_bracket_count_among_four": len(verified),
                "first_verified_bracket_rank": int(verified.proposal_rank.min()) if len(verified) else np.nan,
                "recovered_exact_edge_neighborhood": bool(edge_row.recovered_exact_edge_neighborhood),
                "minimum_exact_edge_orbit_distance": edge_row.minimum_edge_orbit_distance,
                "edge_shadow_duration": edge_row.edge_shadow_duration,
                "direct_period_estimate": edge_row.direct_period_estimate,
                "period_relative_error": edge_row.period_relative_error,
                "lengthscale_1": hyper["lengthscale_1"],
                "lengthscale_2": hyper["lengthscale_2"],
                "signal_std": hyper["signal_std"],
                "hyperparameter_bound_flag": (
                    hyper["signal_std"] >= 0.999999 * 30.0
                    or hyper["lengthscale_1"] >= 0.999999 * 3.0
                    or hyper["lengthscale_2"] >= 0.999999 * 3.0
                ),
            }
        )
    comparison = pd.DataFrame(rows)
    band_ratio = (
        uncertainty["active"]["frozen_separator_band_normalized_integral_reduction_per_added_trajectory"]
        / uncertainty["baseline"]["frozen_separator_band_normalized_integral_reduction_per_added_trajectory"]
    )
    summary = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "frozen_model_manifest_hash": sha256_file(MODEL_DIR / "frozen_model_manifest.json"),
        "models_frozen_before_test": True,
        "shared_seed_count": 25,
        "added_budget_per_model": 8,
        "unique_added_PDE_runs": 15,
        "fixed_test_count": 16,
        "fixed_test_B_count": int(np.sum(cases.true_label == "B")),
        "fixed_test_C_count": int(np.sum(cases.true_label == "C")),
        "both_models_hard_label_accuracy": "16/16; not discriminating at this sample size",
        "active_vs_baseline_Brier_ratio": float(
            metrics["active"]["Brier_score_on_BC_cases"]
            / metrics["space_filling_baseline"]["Brier_score_on_BC_cases"]
        ),
        "active_vs_baseline_log_loss_ratio": float(
            metrics["active"]["log_loss_on_BC_cases"]
            / metrics["space_filling_baseline"]["log_loss_on_BC_cases"]
        ),
        "active_vs_baseline_separator_band_uncertainty_reduction_per_run_ratio": float(band_ratio),
        "active_verified_brackets_among_four": int(
            comparison.loc[comparison.model == "active", "verified_bracket_count_among_four"].iloc[0]
        ),
        "baseline_verified_brackets_among_four": int(
            comparison.loc[comparison.model == "baseline", "verified_bracket_count_among_four"].iloc[0]
        ),
        "baseline_edge_recovery": {
            "minimum_orbit_distance": float(edge.loc[edge.model == "baseline", "minimum_edge_orbit_distance"].iloc[0]),
            "shadow_duration_below_0.02": float(edge.loc[edge.model == "baseline", "edge_shadow_duration"].iloc[0]),
            "direct_period_estimate": float(edge.loc[edge.model == "baseline", "direct_period_estimate"].iloc[0]),
            "exact_period": float(edge.loc[edge.model == "baseline", "exact_edge_period"].iloc[0]),
            "period_relative_error": float(edge.loc[edge.model == "baseline", "period_relative_error"].iloc[0]),
        },
        "interpretation": [
            "Active sampling improved proper probabilistic scores and reduced the predeclared uncertainty functional more efficiently on this plane.",
            "Those gains did not translate into a verified B-C bracket among four frozen active proposals; the space-filling rank-4 bracket recovered E_BC.",
            "This single small experiment does not establish that space-filling is generally better for edge discovery; it shows that classification and posterior uncertainty alone are insufficient dynamical-usefulness criteria.",
            "The active signal standard deviation and baseline first length scale reached predeclared optimization bounds, so calibration and hyperparameter-prior sensitivity remain required robustness checks.",
        ],
        "fixed_test_labels_used_for_refitting": False,
        "exact_edge_used_for_fitting_or_acquisition": False,
        "frozen_artifact_hashes_verified": True,
    }

    OUTDIR.mkdir(parents=True)
    comparison.to_csv(OUTDIR / "controlled_model_comparison.csv", index=False)
    brackets.to_csv(OUTDIR / "frozen_bracket_outcomes.csv", index=False)
    with (OUTDIR / "controlled_experiment_summary.json").open("w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2)
    make_figure(OUTDIR / "controlled_experiment_summary.png", cases, metrics, uncertainty, brackets)
    print(comparison.to_string(index=False))
    print(json.dumps(summary, indent=2))
    print(f"Consolidated controlled experiment saved in: {OUTDIR}")


if __name__ == "__main__":
    main()
