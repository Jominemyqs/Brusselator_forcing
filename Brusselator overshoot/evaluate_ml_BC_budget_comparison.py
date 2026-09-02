#!/usr/bin/env python3
"""Open the fixed test once and evaluate the two already-frozen GP models."""

import hashlib
import json
import os
import shutil
from datetime import datetime, timezone
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/private/tmp/brusselator_matplotlib")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist


ROOT = Path(__file__).resolve().parent
MODEL_DIR = ROOT / "experiment_outputs" / "ml_BC_budget_comparison_v1"
TEST_DIR = ROOT / "experiment_outputs" / "ml_BC_fixed_test_labeling_v1"
OUTDIR = ROOT / "experiment_outputs" / "ml_BC_budget_evaluation_v1"


def sha256_file(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def verify_frozen_files(manifest):
    paths = {
        "active_model_hash": MODEL_DIR / "active_probit_gp_model.npz",
        "baseline_model_hash": MODEL_DIR / "baseline_probit_gp_model.npz",
        "active_blind_prediction_hash": MODEL_DIR / "active_fixed_test_blind_predictions.csv",
        "baseline_blind_prediction_hash": MODEL_DIR / "baseline_fixed_test_blind_predictions.csv",
    }
    failed = [key for key, path in paths.items() if sha256_file(path) != manifest[key]]
    if failed:
        raise RuntimeError("A frozen model artifact changed after test access: " + ", ".join(failed))


def model_metrics(true_label, probability):
    binary = np.isin(true_label, ["B", "C"])
    target = (true_label[binary] == "C").astype(float)
    p = np.clip(probability[binary], 1e-15, 1.0 - 1e-15)
    predicted = np.where(probability >= 0.5, "C", "B")
    return {
        "test_count": int(len(true_label)),
        "binary_BC_test_count": int(np.sum(binary)),
        "out_of_taxonomy_count": int(np.sum(~binary)),
        "correct_count": int(np.sum((predicted == true_label) & binary)),
        "accuracy_all_16": float(np.mean((predicted == true_label) & binary)),
        "Brier_score_on_BC_cases": float(np.mean((p - target) ** 2)),
        "log_loss_on_BC_cases": float(-np.mean(target * np.log(p) + (1.0 - target) * np.log(1.0 - p))),
        "mean_probability_C": float(np.mean(probability)),
    }


def distance_to_contour(frame, contour):
    features = ["normalized_alpha_1", "normalized_alpha_2"]
    distances = cdist(frame[features].to_numpy(float), contour[features].to_numpy(float))
    return np.min(distances, axis=1)


def contour_pair_metrics(left, right):
    features = ["normalized_alpha_1", "normalized_alpha_2"]
    distances = cdist(left[features].to_numpy(float), right[features].to_numpy(float))
    return {
        "mean_left_to_right": float(np.mean(np.min(distances, axis=1))),
        "mean_right_to_left": float(np.mean(np.min(distances, axis=0))),
        "Hausdorff_distance": float(max(np.max(np.min(distances, axis=1)), np.max(np.min(distances, axis=0)))),
    }


def make_prediction_figure(path, cases, active_contour, baseline_contour):
    figure, axes = plt.subplots(1, 2, figsize=(12.5, 4.8), constrained_layout=True)
    for axis, prefix, contour, title in [
        (axes[0], "active", active_contour, "Active acquisition"),
        (axes[1], "baseline", baseline_contour, "Space-filling control"),
    ]:
        for _, segment in contour.groupby("segment_id"):
            axis.plot(segment.alpha_1, segment.alpha_2, color="black", linewidth=1.8)
        for label, color, marker in [("B", "tab:blue", "o"), ("C", "tab:orange", "o"), ("other", "tab:purple", "X"), ("U", "0.4", "X")]:
            selected = cases.true_label == label
            if selected.any():
                edge = np.where(cases.loc[selected, f"{prefix}_correct"], "black", "red")
                axis.scatter(
                    cases.loc[selected, "alpha_1"], cases.loc[selected, "alpha_2"],
                    c=color, marker=marker, s=70, edgecolors=edge, linewidths=1.3,
                    label=label, zorder=4
                )
        axis.set_title(title + " (red edge = error)")
        axis.set_xlabel(r"$\alpha_1$")
        axis.set_ylabel(r"$\alpha_2$")
        axis.set_aspect("equal")
        axis.legend(loc="best")
    figure.suptitle("Blind 16-state test against frozen $p_C=0.5$ separators")
    figure.savefig(path, dpi=300)
    plt.close(figure)


def make_uncertainty_figure(path, uncertainty):
    models = ["initial", "active", "baseline"]
    full = [uncertainty[name]["full_plane_normalized_integral"] for name in models]
    band = [uncertainty[name]["frozen_separator_band_normalized_integral"] for name in models]
    x = np.arange(3)
    figure, axes = plt.subplots(1, 2, figsize=(9.5, 4.2), constrained_layout=True)
    axes[0].bar(x, full, color=["0.5", "tab:green", "tab:cyan"])
    axes[1].bar(x, band, color=["0.5", "tab:green", "tab:cyan"])
    for axis, title in zip(axes, ["Full frozen plane", "Frozen initial-separator band"]):
        axis.set_xticks(x, models)
        axis.set_ylabel(r"$\int p_C(1-p_C)\,dz$")
        axis.set_title(title)
    figure.suptitle("Posterior uncertainty after equal eight-trajectory budgets")
    figure.savefig(path, dpi=300)
    plt.close(figure)


def main():
    if OUTDIR.exists():
        raise FileExistsError(f"Refusing to overwrite existing evaluation: {OUTDIR}")
    frozen = json.loads((MODEL_DIR / "frozen_model_manifest.json").read_text(encoding="utf-8"))
    if not frozen["models_frozen"] or frozen["fixed_test_labels_used"]:
        raise RuntimeError("The comparison was not frozen in the required blind state")
    verify_frozen_files(frozen)

    labels = pd.read_csv(TEST_DIR / "fixed_test_labels.csv")
    active = pd.read_csv(MODEL_DIR / "active_fixed_test_blind_predictions.csv")
    baseline = pd.read_csv(MODEL_DIR / "baseline_fixed_test_blind_predictions.csv")
    if list(labels.state_index) != frozen["fixed_test_state_indices"]:
        raise RuntimeError("Fixed-test state membership or order changed after model freezing")
    if set(active.columns) & {"final_label", "outcome_label", "true_label"}:
        raise RuntimeError("The active blind prediction file contains an outcome label")
    if set(baseline.columns) & {"final_label", "outcome_label", "true_label"}:
        raise RuntimeError("The baseline blind prediction file contains an outcome label")

    cases = labels[[
        "state_index", "sample_id", "alpha_1", "alpha_2",
        "normalized_alpha_1", "normalized_alpha_2", "final_label",
        "raw_five_way_outcome", "tested_duration"
    ]].rename(columns={"final_label": "true_label"})
    cases = cases.merge(
        active[["state_index", "probability_C", "predictive_entropy_bits"]].rename(
            columns={"probability_C": "active_probability_C", "predictive_entropy_bits": "active_entropy_bits"}
        ), on="state_index", validate="one_to_one"
    ).merge(
        baseline[["state_index", "probability_C", "predictive_entropy_bits"]].rename(
            columns={"probability_C": "baseline_probability_C", "predictive_entropy_bits": "baseline_entropy_bits"}
        ), on="state_index", validate="one_to_one"
    )
    for name in ["active", "baseline"]:
        cases[f"{name}_predicted_label"] = np.where(cases[f"{name}_probability_C"] >= 0.5, "C", "B")
        cases[f"{name}_correct"] = cases[f"{name}_predicted_label"] == cases.true_label

    active_contour = pd.read_csv(MODEL_DIR / "active_probability_half_contour.csv")
    baseline_contour = pd.read_csv(MODEL_DIR / "baseline_probability_half_contour.csv")
    initial_contour = pd.read_csv(MODEL_DIR / "initial_probability_half_contour.csv")
    cases["active_distance_to_separator_normalized"] = distance_to_contour(cases, active_contour)
    cases["baseline_distance_to_separator_normalized"] = distance_to_contour(cases, baseline_contour)

    active_metrics = model_metrics(cases.true_label.to_numpy(), cases.active_probability_C.to_numpy(float))
    baseline_metrics = model_metrics(cases.true_label.to_numpy(), cases.baseline_probability_C.to_numpy(float))
    both_correct = int(np.sum(cases.active_correct & cases.baseline_correct))
    active_only = int(np.sum(cases.active_correct & ~cases.baseline_correct))
    baseline_only = int(np.sum(~cases.active_correct & cases.baseline_correct))
    both_wrong = int(np.sum(~cases.active_correct & ~cases.baseline_correct))
    uncertainty = json.loads((MODEL_DIR / "uncertainty_metrics.json").read_text(encoding="utf-8"))
    contour_metrics = {
        "initial_vs_active": contour_pair_metrics(initial_contour, active_contour),
        "initial_vs_baseline": contour_pair_metrics(initial_contour, baseline_contour),
        "active_vs_baseline": contour_pair_metrics(active_contour, baseline_contour),
    }
    summary = {
        "evaluation_created_utc": datetime.now(timezone.utc).isoformat(),
        "models_were_frozen_before_test_labels": True,
        "fixed_test_count": len(cases),
        "active": active_metrics,
        "space_filling_baseline": baseline_metrics,
        "paired_correctness": {
            "both_correct": both_correct,
            "active_only_correct": active_only,
            "baseline_only_correct": baseline_only,
            "both_wrong": both_wrong,
            "interpretation": "with n=16, one-outcome differences are descriptive rather than strong statistical evidence",
        },
        "separator_comparison_normalized_coordinates": contour_metrics,
        "uncertainty_functional": uncertainty,
        "dynamical_usefulness_status": "not assessed by classification metrics; ML bracket proposals require full-PDE edge tracking",
    }

    OUTDIR.mkdir(parents=True)
    cases.to_csv(OUTDIR / "blind_test_individual_predictions.csv", index=False)
    shutil.copy2(MODEL_DIR / "frozen_model_manifest.json", OUTDIR / "evaluated_frozen_model_manifest.json")
    with (OUTDIR / "blind_test_metrics.json").open("w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2)
    make_prediction_figure(OUTDIR / "blind_test_separator_predictions.png", cases, active_contour, baseline_contour)
    make_uncertainty_figure(OUTDIR / "budget_matched_uncertainty.png", uncertainty)
    print(json.dumps(summary, indent=2))
    print(f"Blind evaluation saved in: {OUTDIR}")


if __name__ == "__main__":
    main()
