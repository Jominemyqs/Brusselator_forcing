#!/usr/bin/env python3
"""Refit and freeze budget-matched B--C probit-GP models.

This program intentionally has no input path for fixed-test outcomes.  It writes
the two models and their blind predictions before the test-label runner is
permitted to execute.
"""

import hashlib
import json
import os
import shutil
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/private/tmp/brusselator_matplotlib")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist

from probit_gp import optimize_hyperparameters, predict_proba


ROOT = Path(__file__).resolve().parent
CONFIG_FILE = ROOT / "ml_bc_budget_comparison_config.json"


def sha256_file(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def predictive_entropy(probability):
    p = np.clip(np.asarray(probability, dtype=float), 1e-15, 1.0 - 1e-15)
    return -(p * np.log2(p) + (1.0 - p) * np.log2(1.0 - p))


def fit_model(training, features, config):
    X = training[features].to_numpy(float)
    y = np.where(training["final_label"].to_numpy() == config["positive_class"], 1.0, -1.0)
    kernel = config["kernel"]
    laplace = config["laplace"]
    optimizer = config["optimizer"]
    state, records = optimize_hyperparameters(
        X,
        y,
        kernel["initial_starts"],
        lengthscale_bounds=tuple(kernel["lengthscale_bounds"]),
        signal_std_bounds=tuple(kernel["signal_std_bounds"]),
        jitter=float(kernel["jitter"]),
        maximum_laplace_iterations=int(laplace["maximum_iterations"]),
        laplace_tolerance=float(laplace["tolerance"]),
        maximum_optimizer_iterations=int(optimizer["maximum_iterations"]),
        function_tolerance=float(optimizer["function_tolerance"]),
    )
    if not state.converged:
        raise RuntimeError("The Laplace iteration did not converge")
    return state, pd.DataFrame(records)


def add_predictions(frame, state, features):
    result = frame.copy()
    probability, mean, variance = predict_proba(state, result[features].to_numpy(float))
    result["probability_C"] = probability
    result["probability_B"] = 1.0 - probability
    result["predictive_entropy_bits"] = predictive_entropy(probability)
    result["posterior_variance_functional"] = probability * (1.0 - probability)
    result["latent_mean"] = mean
    result["latent_variance"] = variance
    return result


def contour_frame(axis, probability):
    figure, plotting_axis = plt.subplots()
    contour = plotting_axis.contour(axis, axis, probability, levels=[0.5])
    rows = []
    for segment_id, segment in enumerate(contour.allsegs[0], start=1):
        for point_id, point in enumerate(segment, start=1):
            rows.append(
                {
                    "segment_id": segment_id,
                    "point_index": point_id,
                    "normalized_alpha_1": float(point[0]),
                    "normalized_alpha_2": float(point[1]),
                }
            )
    plt.close(figure)
    return pd.DataFrame(rows)


def contour_distance(left, right):
    if left.empty or right.empty:
        return {"mean_left_to_right": np.nan, "mean_right_to_left": np.nan, "hausdorff": np.nan}
    columns = ["normalized_alpha_1", "normalized_alpha_2"]
    distances = cdist(left[columns].to_numpy(float), right[columns].to_numpy(float))
    return {
        "mean_left_to_right": float(np.mean(np.min(distances, axis=1))),
        "mean_right_to_left": float(np.mean(np.min(distances, axis=0))),
        "hausdorff": float(max(np.max(np.min(distances, axis=1)), np.max(np.min(distances, axis=0)))),
    }


def uncertainty_metrics(probability, initial_probability, spacing, alpha_scales, margin):
    uncertainty = probability * (1.0 - probability)
    frozen_mask = np.abs(initial_probability - 0.5) <= margin
    cell_area = spacing**2
    one_dimensional_weights = np.ones(probability.shape[0])
    one_dimensional_weights[[0, -1]] = 0.5
    quadrature_weights = np.outer(one_dimensional_weights, one_dimensional_weights)
    physical_factor = float(alpha_scales[0] * alpha_scales[1])
    return {
        "full_plane_normalized_integral": float(np.sum(uncertainty * quadrature_weights) * cell_area),
        "full_plane_physical_integral": float(np.sum(uncertainty * quadrature_weights) * cell_area * physical_factor),
        "frozen_separator_band_normalized_integral": float(np.sum(uncertainty * quadrature_weights * frozen_mask) * cell_area),
        "frozen_separator_band_physical_integral": float(np.sum(uncertainty * quadrature_weights * frozen_mask) * cell_area * physical_factor),
        "frozen_separator_band_mean": float(
            np.sum(uncertainty * quadrature_weights * frozen_mask)
            / np.sum(quadrature_weights * frozen_mask)
        ),
        "frozen_separator_band_grid_count": int(np.sum(frozen_mask)),
    }


def edge_bracket_proposals(manifest, state, features, excluded_indices, proposal_count):
    """Select label-blind adjacent pairs predicted on opposite sides of p_C=0.5."""
    eligible = manifest[
        (manifest.design_role == "acquisition_pool")
        & ~manifest.state_index.isin(excluded_indices)
    ].copy()
    eligible = add_predictions(eligible, state, features)
    lookup = {
        (int(row.row_index), int(row.column_index)): row
        for row in eligible.itertuples(index=False)
    }
    pairs = []
    for (row_index, column_index), left in lookup.items():
        for neighbor_key in [(row_index + 1, column_index), (row_index, column_index + 1)]:
            if neighbor_key not in lookup:
                continue
            right = lookup[neighbor_key]
            left_side = left.probability_C >= 0.5
            right_side = right.probability_C >= 0.5
            if left_side == right_side:
                continue
            B = left if not left_side else right
            C = right if not left_side else left
            separation = float(
                np.hypot(
                    B.normalized_alpha_1 - C.normalized_alpha_1,
                    B.normalized_alpha_2 - C.normalized_alpha_2,
                )
            )
            pairs.append(
                {
                    "B_state_index": int(B.state_index),
                    "C_state_index": int(C.state_index),
                    "B_sample_id": B.sample_id,
                    "C_sample_id": C.sample_id,
                    "B_probability_C": float(B.probability_C),
                    "C_probability_C": float(C.probability_C),
                    "maximum_probability_margin": float(
                        max(abs(B.probability_C - 0.5), abs(C.probability_C - 0.5))
                    ),
                    "sum_probability_margin": float(
                        abs(B.probability_C - 0.5) + abs(C.probability_C - 0.5)
                    ),
                    "normalized_endpoint_separation": separation,
                }
            )
    candidates = pd.DataFrame(pairs).sort_values(
        ["maximum_probability_margin", "sum_probability_margin", "normalized_endpoint_separation", "B_state_index", "C_state_index"]
    )
    selected = []
    used = set()
    for row in candidates.itertuples(index=False):
        endpoints = {int(row.B_state_index), int(row.C_state_index)}
        if endpoints & used:
            continue
        selected.append(row._asdict())
        used |= endpoints
        if len(selected) == proposal_count:
            break
    if len(selected) != proposal_count:
        raise RuntimeError("Could not find the requested number of disjoint predicted edge brackets")
    result = pd.DataFrame(selected)
    result.insert(0, "proposal_rank", np.arange(1, len(result) + 1))
    return result


def save_model(path, state, training):
    np.savez_compressed(
        path,
        X=state.X,
        y=state.y,
        state_indices=training["state_index"].to_numpy(int),
        lengthscales=state.lengthscales,
        signal_std=np.array([state.signal_std]),
        jitter=np.array([state.jitter]),
        f=state.f,
        a=state.a,
        sqrt_W=state.sqrt_W,
        L=state.L,
    )


def git_commit():
    try:
        result = subprocess.run(
            ["git", "rev-parse", "HEAD"], cwd=ROOT, check=True, capture_output=True, text=True
        )
        return result.stdout.strip()
    except (subprocess.CalledProcessError, FileNotFoundError):
        return ""


def make_figure(path, manifest, training_sets, probability_grids, contours, axis, alpha_scales):
    physical_axis_1 = axis * alpha_scales[0]
    physical_axis_2 = axis * alpha_scales[1]
    figure, axes = plt.subplots(1, 3, figsize=(15.5, 4.7), constrained_layout=True)
    titles = {
        "initial": "Initial: 25 labels",
        "active": "Active: 25 + 8 labels",
        "baseline": "Space-filling: 25 + 8 labels",
    }
    for panel, key in zip(axes, ["initial", "active", "baseline"]):
        image = panel.pcolormesh(
            physical_axis_1,
            physical_axis_2,
            probability_grids[key],
            shading="auto",
            vmin=0,
            vmax=1,
            cmap="coolwarm",
        )
        for contour_key, color, linestyle in [
            ("initial", "black", "--"),
            (key, "white", "-"),
        ]:
            frame = contours[contour_key]
            for _, segment in frame.groupby("segment_id"):
                panel.plot(
                    segment["normalized_alpha_1"] * alpha_scales[0],
                    segment["normalized_alpha_2"] * alpha_scales[1],
                    color=color,
                    linestyle=linestyle,
                    linewidth=1.6,
                )
        training = training_sets[key]
        if training is not None:
            for label, color in [("B", "tab:blue"), ("C", "tab:orange")]:
                selected = training["final_label"] == label
                panel.scatter(
                    training.loc[selected, "alpha_1"],
                    training.loc[selected, "alpha_2"],
                    s=34,
                    c=color,
                    edgecolors="black",
                    linewidths=0.35,
                    zorder=4,
                )
        fixed = manifest[manifest["design_role"] == "fixed_test"]
        panel.scatter(
            fixed["alpha_1"], fixed["alpha_2"], marker="s", s=34,
            facecolors="none", edgecolors="black", linewidths=0.7, zorder=3
        )
        panel.set_title(titles[key])
        panel.set_xlabel(r"$\alpha_1$")
        panel.set_ylabel(r"$\alpha_2$")
        panel.set_aspect("equal")
    colorbar = figure.colorbar(image, ax=axes, shrink=0.82)
    colorbar.set_label(r"$p_C$")
    figure.suptitle("Budget-matched models; fixed-test squares remain outcome-blind")
    figure.savefig(path, dpi=300)
    plt.close(figure)


def main():
    started = time.perf_counter()
    with CONFIG_FILE.open("r", encoding="utf-8") as handle:
        config = json.load(handle)
    outdir = ROOT / config["output_root"] / config["experiment_name"]
    if outdir.exists():
        raise FileExistsError(f"Refusing to overwrite existing output: {outdir}")

    manifest = pd.read_csv(ROOT / config["geometry_manifest_csv"])
    seed_labels = pd.read_csv(ROOT / config["initial_labels_csv"])
    acquired = pd.read_csv(ROOT / config["acquisition_labels_csv"])
    active_proposals = pd.read_csv(ROOT / config["active_proposals_csv"])
    baseline_proposals = pd.read_csv(ROOT / config["baseline_proposals_csv"])
    features = config["feature_columns"]

    seed = manifest.merge(
        seed_labels[["state_index", "final_label"]], on="state_index", validate="one_to_one"
    ).sort_values("state_index")
    if len(seed) != config["seed_budget"] or set(seed["final_label"]) - {"B", "C"}:
        raise ValueError("The shared seed set must contain exactly 25 resolved B/C labels")
    if not np.all(seed["design_role"] == "initial_seed"):
        raise ValueError("A non-seed state leaked into the shared training set")
    if len(active_proposals) != config["added_budget_per_model"] or len(baseline_proposals) != config["added_budget_per_model"]:
        raise ValueError("The active and baseline additions must each have the same eight-state budget")
    if len(set(active_proposals.state_index) & set(baseline_proposals.state_index)) != 1:
        raise ValueError("The frozen proposal overlap has changed")
    if set(acquired["final_label"]) - {"B", "C"}:
        raise ValueError("A selected acquisition state did not resolve to B or C; do not substitute another state")

    acquired_core = manifest.merge(
        acquired[["state_index", "final_label"]], on="state_index", validate="one_to_one"
    )
    active_added = acquired_core[acquired_core.state_index.isin(active_proposals.state_index)].copy()
    baseline_added = acquired_core[acquired_core.state_index.isin(baseline_proposals.state_index)].copy()
    if len(active_added) != 8 or len(baseline_added) != 8:
        raise ValueError("Not every frozen acquisition proposal has a completed label")
    active_training = pd.concat([seed, active_added], ignore_index=True).sort_values("state_index")
    baseline_training = pd.concat([seed, baseline_added], ignore_index=True).sort_values("state_index")
    if len(active_training) != 33 or len(baseline_training) != 33:
        raise ValueError("Each comparison model must use exactly 33 trajectories")

    active_state, active_optimization = fit_model(active_training, features, config)
    baseline_state, baseline_optimization = fit_model(baseline_training, features, config)

    fixed_test = manifest[manifest["design_role"] == "fixed_test"].copy()
    if len(fixed_test) != config["fixed_test_count"]:
        raise ValueError("The frozen test set no longer has 16 states")
    if fixed_test["outcome_label"].notna().any():
        raise ValueError("A fixed-test label leaked into the geometry manifest")
    blind_columns = [
        "state_index", "sample_id", "alpha_1", "alpha_2",
        "normalized_alpha_1", "normalized_alpha_2", "design_role"
    ]
    active_blind = add_predictions(fixed_test[blind_columns], active_state, features)
    baseline_blind = add_predictions(fixed_test[blind_columns], baseline_state, features)

    grid_size = int(config["dense_grid_size"])
    axis = np.linspace(-1.0, 1.0, grid_size)
    grid_1, grid_2 = np.meshgrid(axis, axis)
    dense_X = np.column_stack([grid_1.ravel(), grid_2.ravel()])
    initial_dense = pd.read_csv(ROOT / config["initial_model_directory"] / "dense_probability_grid.csv")
    initial_probability = initial_dense["probability_C"].to_numpy(float).reshape(grid_1.shape)
    active_probability = predict_proba(active_state, dense_X)[0].reshape(grid_1.shape)
    baseline_probability = predict_proba(baseline_state, dense_X)[0].reshape(grid_1.shape)
    probability_grids = {
        "initial": initial_probability,
        "active": active_probability,
        "baseline": baseline_probability,
    }
    alpha_scales = [float(manifest.alpha_1.abs().max()), float(manifest.alpha_2.abs().max())]
    contours = {
        "initial": pd.read_csv(ROOT / config["initial_model_directory"] / "probability_half_contour.csv"),
        "active": contour_frame(axis, active_probability),
        "baseline": contour_frame(axis, baseline_probability),
    }
    for frame in contours.values():
        frame["alpha_1"] = frame["normalized_alpha_1"] * alpha_scales[0]
        frame["alpha_2"] = frame["normalized_alpha_2"] * alpha_scales[1]

    spacing = float(axis[1] - axis[0])
    margin = float(config["frozen_separator_band"]["maximum_probability_margin"])
    uncertainty = {
        name: uncertainty_metrics(probability, initial_probability, spacing, alpha_scales, margin)
        for name, probability in probability_grids.items()
    }
    for model_name in ["active", "baseline"]:
        for domain in ["full_plane_normalized_integral", "frozen_separator_band_normalized_integral"]:
            initial_value = uncertainty["initial"][domain]
            final_value = uncertainty[model_name][domain]
            uncertainty[model_name][domain + "_reduction"] = float(initial_value - final_value)
            uncertainty[model_name][domain + "_reduction_per_added_trajectory"] = float(
                (initial_value - final_value) / config["added_budget_per_model"]
            )

    contour_comparison = {}
    for left, right in [("initial", "active"), ("initial", "baseline"), ("active", "baseline")]:
        contour_comparison[f"{left}_to_{right}"] = contour_distance(contours[left], contours[right])

    target_285 = manifest[manifest.state_index == config["state_285_index"]]
    if len(target_285) != 1:
        raise ValueError("The frozen state-285 diagnostic is missing")
    state_285 = {
        "state_index": int(config["state_285_index"]),
        "alpha_1": float(target_285.iloc[0].alpha_1),
        "alpha_2": float(target_285.iloc[0].alpha_2),
    }
    initial_training_predictions = pd.read_csv(
        ROOT / config["initial_model_directory"] / "training_posterior_predictions.csv"
    )
    state_285["initial_probability_C"] = float(
        initial_training_predictions.loc[
            initial_training_predictions.state_index == config["state_285_index"],
            "probability_C",
        ].iloc[0]
    )
    initial_loo = pd.read_csv(
        ROOT / config["initial_model_directory"] / "nested_LOO_predictions.csv"
    )
    state_285["initial_leave_one_out_probability_C"] = float(
        initial_loo.loc[
            initial_loo.state_index == config["state_285_index"], "probability_C"
        ].iloc[0]
    )
    state_285["active_probability_C"] = float(predict_proba(active_state, target_285[features].to_numpy(float))[0][0])
    state_285["baseline_probability_C"] = float(predict_proba(baseline_state, target_285[features].to_numpy(float))[0][0])
    point_285 = target_285[features].to_numpy(float)
    for model_name, contour in contours.items():
        contour_points = contour[features].to_numpy(float)
        state_285[f"distance_to_{model_name}_separator_normalized"] = float(
            np.min(cdist(point_285, contour_points))
        )

    excluded = set(seed.state_index) | set(acquired_core.state_index) | set(fixed_test.state_index)
    bracket_count = int(config["edge_bracket_proposals_per_model"])
    active_brackets = edge_bracket_proposals(
        manifest, active_state, features, excluded, bracket_count
    )
    baseline_brackets = edge_bracket_proposals(
        manifest, baseline_state, features, excluded, bracket_count
    )

    outdir.mkdir(parents=True)
    shutil.copy2(CONFIG_FILE, outdir / CONFIG_FILE.name)
    active_training.to_csv(outdir / "active_training_set.csv", index=False)
    baseline_training.to_csv(outdir / "baseline_training_set.csv", index=False)
    active_optimization.to_csv(outdir / "active_hyperparameter_optimization.csv", index=False)
    baseline_optimization.to_csv(outdir / "baseline_hyperparameter_optimization.csv", index=False)
    active_blind.to_csv(outdir / "active_fixed_test_blind_predictions.csv", index=False)
    baseline_blind.to_csv(outdir / "baseline_fixed_test_blind_predictions.csv", index=False)
    active_brackets.to_csv(outdir / "active_edge_bracket_proposals.csv", index=False)
    baseline_brackets.to_csv(outdir / "baseline_edge_bracket_proposals.csv", index=False)
    for name, probability in probability_grids.items():
        dense = pd.DataFrame(
            {
                "normalized_alpha_1": dense_X[:, 0],
                "normalized_alpha_2": dense_X[:, 1],
                "alpha_1": dense_X[:, 0] * alpha_scales[0],
                "alpha_2": dense_X[:, 1] * alpha_scales[1],
                "probability_C": probability.ravel(),
                "probability_B": 1.0 - probability.ravel(),
                "posterior_variance_functional": (probability * (1.0 - probability)).ravel(),
            }
        )
        dense.to_csv(outdir / f"{name}_dense_probability_grid.csv", index=False)
        contours[name].to_csv(outdir / f"{name}_probability_half_contour.csv", index=False)
    with (outdir / "uncertainty_metrics.json").open("w", encoding="utf-8") as handle:
        json.dump(uncertainty, handle, indent=2)
    with (outdir / "contour_comparison.json").open("w", encoding="utf-8") as handle:
        json.dump(contour_comparison, handle, indent=2)
    with (outdir / "state_285_diagnostic.json").open("w", encoding="utf-8") as handle:
        json.dump(state_285, handle, indent=2)

    active_model_path = outdir / "active_probit_gp_model.npz"
    baseline_model_path = outdir / "baseline_probit_gp_model.npz"
    save_model(active_model_path, active_state, active_training)
    save_model(baseline_model_path, baseline_state, baseline_training)
    make_figure(
        outdir / "budget_matched_separator_comparison.png",
        manifest,
        {"initial": seed, "active": active_training, "baseline": baseline_training},
        probability_grids,
        contours,
        axis,
        alpha_scales,
    )

    summary = {
        "experiment_name": config["experiment_name"],
        "models_frozen": True,
        "fixed_test_labels_used": False,
        "shared_seed_count": len(seed),
        "active_added_count": len(active_added),
        "baseline_added_count": len(baseline_added),
        "active_training_count": len(active_training),
        "baseline_training_count": len(baseline_training),
        "proposal_overlap_count": len(set(active_added.state_index) & set(baseline_added.state_index)),
        "active_hyperparameters": {
            "lengthscale_1": float(active_state.lengthscales[0]),
            "lengthscale_2": float(active_state.lengthscales[1]),
            "signal_std": float(active_state.signal_std),
        },
        "baseline_hyperparameters": {
            "lengthscale_1": float(baseline_state.lengthscales[0]),
            "lengthscale_2": float(baseline_state.lengthscales[1]),
            "signal_std": float(baseline_state.signal_std),
        },
        "uncertainty": uncertainty,
        "contour_comparison": contour_comparison,
        "state_285": state_285,
        "active_edge_bracket_proposal_count": len(active_brackets),
        "baseline_edge_bracket_proposal_count": len(baseline_brackets),
        "edge_bracket_endpoints_labeled_during_selection": False,
        "exact_edge_orbit_used_for_fitting": False,
        "held_out_CA_region_used": False,
    }
    with (outdir / "model_summary.json").open("w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2)

    frozen_manifest = {
        "schema_version": config["schema_version"],
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "models_frozen": True,
        "fixed_test_labels_used": False,
        "active_model_hash": sha256_file(active_model_path),
        "baseline_model_hash": sha256_file(baseline_model_path),
        "active_blind_prediction_hash": sha256_file(outdir / "active_fixed_test_blind_predictions.csv"),
        "baseline_blind_prediction_hash": sha256_file(outdir / "baseline_fixed_test_blind_predictions.csv"),
        "active_edge_bracket_proposal_hash": sha256_file(outdir / "active_edge_bracket_proposals.csv"),
        "baseline_edge_bracket_proposal_hash": sha256_file(outdir / "baseline_edge_bracket_proposals.csv"),
        "configuration_hash": sha256_file(outdir / CONFIG_FILE.name),
        "fixed_test_state_indices": [int(value) for value in active_blind.state_index],
    }
    with (outdir / "frozen_model_manifest.json").open("w", encoding="utf-8") as handle:
        json.dump(frozen_manifest, handle, indent=2)

    metadata = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "runtime_seconds": time.perf_counter() - started,
        "python_version": sys.version,
        "numpy_version": np.__version__,
        "pandas_version": pd.__version__,
        "matplotlib_version": matplotlib.__version__,
        "git_commit": git_commit(),
    }
    with (outdir / "run_metadata.json").open("w", encoding="utf-8") as handle:
        json.dump(metadata, handle, indent=2)
    print(json.dumps(summary, indent=2))
    print(f"Frozen budget-matched models saved in: {outdir}")


if __name__ == "__main__":
    main()
