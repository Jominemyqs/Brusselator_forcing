#!/usr/bin/env python3
"""Fit the first probabilistic B--C boundary model and freeze proposals."""

import json
import math
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

from probit_gp import optimize_hyperparameters, predict_proba


ROOT = Path(__file__).resolve().parent
CONFIG_FILE = ROOT / "ml_bc_probabilistic_config.json"


def load_config():
    with CONFIG_FILE.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def predictive_entropy(probability):
    p = np.clip(np.asarray(probability, dtype=float), 1e-15, 1.0 - 1e-15)
    return -(p * np.log2(p) + (1.0 - p) * np.log2(1.0 - p))


def minimum_distances(X, reference):
    difference = X[:, None, :] - reference[None, :, :]
    return np.sqrt(np.min(np.sum(difference**2, axis=2), axis=1))


def add_predictions(frame, state, features):
    X = frame[features].to_numpy(float)
    probability, latent_mean, latent_variance = predict_proba(state, X)
    result = frame.copy()
    result["probability_C"] = probability
    result["probability_B"] = 1.0 - probability
    result["predictive_entropy_bits"] = predictive_entropy(probability)
    result["probability_margin"] = np.abs(probability - 0.5)
    result["latent_mean"] = latent_mean
    result["latent_variance"] = latent_variance
    return result


def active_batch(candidate, training_X, config):
    X = candidate[["normalized_alpha_1", "normalized_alpha_2"]].to_numpy(float)
    state_indices = candidate["state_index"].to_numpy(int)
    entropy = candidate["predictive_entropy_bits"].to_numpy(float)
    options = config["active_acquisition"]
    floor = float(options["floor_weight"])
    novelty = minimum_distances(X, training_X)
    novelty_factor = floor + (1.0 - floor) * np.minimum(
        novelty / float(options["novelty_scale"]), 1.0
    )
    available = np.ones(len(candidate), dtype=bool)
    chosen = []
    selection_records = []
    for step in range(1, int(config["active_batch_size"]) + 1):
        if chosen:
            selected_distance = minimum_distances(X, X[np.asarray(chosen)])
            diversity_factor = floor + (1.0 - floor) * np.minimum(
                selected_distance / float(options["diversity_scale"]), 1.0
            )
        else:
            selected_distance = np.full(len(candidate), np.nan)
            diversity_factor = np.ones(len(candidate))
        score = entropy * novelty_factor * diversity_factor
        score[~available] = -np.inf
        order = np.lexsort((state_indices, -score))
        selected = int(order[0])
        chosen.append(selected)
        available[selected] = False
        selection_records.append(
            {
                "state_index": int(state_indices[selected]),
                "acquisition_step": step,
                "acquisition_score": float(score[selected]),
                "training_novelty_distance": float(novelty[selected]),
                "training_novelty_factor": float(novelty_factor[selected]),
                "selected_batch_distance": (
                    math.nan
                    if step == 1
                    else float(selected_distance[selected])
                ),
                "batch_diversity_factor": float(diversity_factor[selected]),
            }
        )
    selection = pd.DataFrame(selection_records)
    return selection.merge(candidate, on="state_index", how="left", validate="one_to_one")


def space_filling_batch(candidate, training_X, batch_size):
    X = candidate[["normalized_alpha_1", "normalized_alpha_2"]].to_numpy(float)
    state_indices = candidate["state_index"].to_numpy(int)
    available = np.ones(len(candidate), dtype=bool)
    reference = training_X.copy()
    records = []
    for step in range(1, int(batch_size) + 1):
        distance = minimum_distances(X, reference)
        distance[~available] = -np.inf
        order = np.lexsort((state_indices, -distance))
        selected = int(order[0])
        records.append(
            {
                "state_index": int(state_indices[selected]),
                "baseline_step": step,
                "space_filling_distance": float(distance[selected]),
            }
        )
        available[selected] = False
        reference = np.vstack([reference, X[selected]])
    selection = pd.DataFrame(records)
    return selection.merge(candidate, on="state_index", how="left", validate="one_to_one")


def fit_model(X, y, config):
    kernel = config["kernel"]
    laplace = config["laplace"]
    optimizer = config["optimizer"]
    return optimize_hyperparameters(
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


def nested_leave_one_out(training, features, config):
    X = training[features].to_numpy(float)
    y = np.where(training["final_label"].to_numpy() == "C", 1.0, -1.0)
    rows = []
    for held_out in range(len(training)):
        keep = np.arange(len(training)) != held_out
        state, _ = fit_model(X[keep], y[keep], config)
        probability, latent_mean, latent_variance = predict_proba(
            state, X[held_out : held_out + 1]
        )
        rows.append(
            {
                "state_index": int(training.iloc[held_out]["state_index"]),
                "true_label": training.iloc[held_out]["final_label"],
                "probability_C": float(probability[0]),
                "predicted_label": "C" if probability[0] >= 0.5 else "B",
                "latent_mean": float(latent_mean[0]),
                "latent_variance": float(latent_variance[0]),
                "lengthscale_1": float(state.lengthscales[0]),
                "lengthscale_2": float(state.lengthscales[1]),
                "signal_std": float(state.signal_std),
                "log_marginal_likelihood": float(state.log_marginal_likelihood),
                "laplace_converged": bool(state.converged),
            }
        )
    return pd.DataFrame(rows)


def metrics_from_loo(loo):
    target = (loo["true_label"].to_numpy() == "C").astype(float)
    probability = np.clip(loo["probability_C"].to_numpy(float), 1e-15, 1 - 1e-15)
    predicted = (probability >= 0.5).astype(float)
    return {
        "nested_LOO_accuracy": float(np.mean(predicted == target)),
        "nested_LOO_Brier_score": float(np.mean((probability - target) ** 2)),
        "nested_LOO_log_loss": float(
            -np.mean(target * np.log(probability) + (1 - target) * np.log(1 - probability))
        ),
    }


def contour_records(axis_1, axis_2, probability_grid):
    figure, axis = plt.subplots()
    contour = axis.contour(axis_1, axis_2, probability_grid, levels=[0.5])
    rows = []
    for segment_id, segment in enumerate(contour.allsegs[0], start=1):
        for point_index, point in enumerate(segment, start=1):
            rows.append(
                {
                    "segment_id": segment_id,
                    "point_index": point_index,
                    "normalized_alpha_1": float(point[0]),
                    "normalized_alpha_2": float(point[1]),
                }
            )
    plt.close(figure)
    return pd.DataFrame(rows)


def make_figure(
    output_file,
    training,
    fixed_test,
    active,
    baseline,
    grid_1,
    grid_2,
    probability,
    entropy,
    alpha_scales,
):
    physical_1 = grid_1 * alpha_scales[0]
    physical_2 = grid_2 * alpha_scales[1]
    figure, axes = plt.subplots(1, 2, figsize=(12.5, 5.0), constrained_layout=True)
    for axis, field, title, color_label in [
        (axes[0], probability, "Posterior probability of C", "$p_C$"),
        (axes[1], entropy, "Predictive entropy", "bits"),
    ]:
        image = axis.pcolormesh(
            physical_1,
            physical_2,
            field,
            shading="auto",
            cmap="viridis",
        )
        axis.contour(physical_1, physical_2, probability, levels=[0.5], colors="white", linewidths=2)
        for label, color, marker in [("B", "tab:blue", "o"), ("C", "tab:orange", "o")]:
            selected = training["final_label"] == label
            axis.scatter(
                training.loc[selected, "alpha_1"],
                training.loc[selected, "alpha_2"],
                s=55,
                c=color,
                marker=marker,
                edgecolors="black",
                linewidths=0.5,
                label=f"training {label}",
                zorder=4,
            )
        axis.scatter(
            fixed_test["alpha_1"],
            fixed_test["alpha_2"],
            s=45,
            facecolors="none",
            edgecolors="white",
            marker="s",
            linewidths=1.2,
            label="blind fixed test",
            zorder=4,
        )
        axis.scatter(
            active["alpha_1"],
            active["alpha_2"],
            s=120,
            c="magenta",
            marker="*",
            edgecolors="black",
            linewidths=0.6,
            label="active proposal",
            zorder=5,
        )
        if axis is axes[1]:
            axis.scatter(
                baseline["alpha_1"],
                baseline["alpha_2"],
                s=75,
                c="cyan",
                marker="x",
                linewidths=1.5,
                label="space-filling baseline",
                zorder=5,
            )
        axis.set_xlabel(r"$\alpha_1$")
        axis.set_ylabel(r"$\alpha_2$")
        axis.set_title(title)
        axis.set_aspect("equal")
        axis.legend(loc="best", fontsize=8)
        colorbar = figure.colorbar(image, ax=axis)
        colorbar.set_label(color_label)
    figure.suptitle(
        "Initial probit-GP boundary model; exact edge and held-out challenge excluded"
    )
    figure.savefig(output_file, dpi=300)
    plt.close(figure)


def git_commit():
    try:
        result = subprocess.run(
            ["git", "rev-parse", "HEAD"],
            cwd=ROOT,
            check=True,
            capture_output=True,
            text=True,
        )
        return result.stdout.strip()
    except (subprocess.CalledProcessError, FileNotFoundError):
        return ""


def main():
    started = time.perf_counter()
    config = load_config()
    outdir = ROOT / config["output_root"] / config["experiment_name"]
    if outdir.exists():
        raise FileExistsError(f"Refusing to overwrite existing output: {outdir}")

    manifest = pd.read_csv(ROOT / config["geometry_manifest_csv"])
    labels = pd.read_csv(ROOT / config["initial_labels_csv"])
    training = manifest.merge(
        labels[["state_index", "final_label"]],
        on="state_index",
        how="inner",
        validate="one_to_one",
    ).sort_values("state_index")
    if len(training) != 25 or set(training["final_label"]) - {"B", "C"}:
        raise ValueError("The first GP requires exactly 25 resolved B/C seed labels")
    if not np.all(training["design_role"] == "initial_seed"):
        raise ValueError("A non-seed state entered the GP training set")
    fixed_test = manifest[manifest["design_role"] == "fixed_test"].copy()
    candidate = manifest[manifest["design_role"] == "acquisition_pool"].copy()
    if len(fixed_test) != 16 or len(candidate) != 248:
        raise ValueError("Frozen test or acquisition pool count has changed")

    features = config["feature_columns"]
    X = training[features].to_numpy(float)
    y = np.where(training["final_label"].to_numpy() == config["positive_class"], 1.0, -1.0)
    state, optimization_records = fit_model(X, y, config)
    if not state.converged:
        raise RuntimeError("Full-data Laplace iteration did not converge")
    training_predictions = add_predictions(training, state, features)
    candidate_predictions = add_predictions(candidate, state, features)
    fixed_test_predictions = add_predictions(fixed_test, state, features)
    fixed_test_predictions = fixed_test_predictions.drop(
        columns=["outcome_label"], errors="ignore"
    )
    candidate_predictions["pure_entropy_rank"] = (
        candidate_predictions["predictive_entropy_bits"]
        .rank(method="first", ascending=False)
        .astype(int)
    )

    active = active_batch(candidate_predictions, X, config)
    baseline = space_filling_batch(candidate_predictions, X, config["active_batch_size"])
    loo = nested_leave_one_out(training, features, config)
    loo_metrics = metrics_from_loo(loo)

    grid_values = np.linspace(-1.0, 1.0, int(config["dense_grid_size"]))
    grid_1, grid_2 = np.meshgrid(grid_values, grid_values)
    dense_X = np.column_stack([grid_1.ravel(), grid_2.ravel()])
    dense_probability, dense_mean, dense_variance = predict_proba(state, dense_X)
    probability_grid = dense_probability.reshape(grid_1.shape)
    entropy_grid = predictive_entropy(probability_grid)
    contour = contour_records(grid_1, grid_2, probability_grid)
    alpha_scales = [float(manifest["alpha_1"].abs().max()), float(manifest["alpha_2"].abs().max())]
    if not contour.empty:
        contour["alpha_1"] = contour["normalized_alpha_1"] * alpha_scales[0]
        contour["alpha_2"] = contour["normalized_alpha_2"] * alpha_scales[1]

    target = (training_predictions["final_label"].to_numpy() == "C").astype(int)
    training_class = (training_predictions["probability_C"].to_numpy() >= 0.5).astype(int)
    summary = {
        "schema_version": config["schema_version"],
        "experiment_name": config["experiment_name"],
        "training_count": int(len(training)),
        "training_B_count": int(np.sum(training["final_label"] == "B")),
        "training_C_count": int(np.sum(training["final_label"] == "C")),
        "fixed_test_count": int(len(fixed_test)),
        "acquisition_pool_count": int(len(candidate)),
        "active_batch_size": int(len(active)),
        "space_filling_baseline_size": int(len(baseline)),
        "lengthscale_1": float(state.lengthscales[0]),
        "lengthscale_2": float(state.lengthscales[1]),
        "signal_std": float(state.signal_std),
        "signal_std_at_upper_bound": bool(
            state.signal_std >= 0.999999 * float(config["kernel"]["signal_std_bounds"][1])
        ),
        "lengthscale_at_bound": bool(
            np.any(state.lengthscales <= 1.000001 * float(config["kernel"]["lengthscale_bounds"][0]))
            or np.any(state.lengthscales >= 0.999999 * float(config["kernel"]["lengthscale_bounds"][1]))
        ),
        "log_marginal_likelihood": float(state.log_marginal_likelihood),
        "laplace_iterations": int(state.iterations),
        "laplace_converged": bool(state.converged),
        "training_accuracy": float(np.mean(training_class == target)),
        **loo_metrics,
        "maximum_candidate_entropy_bits": float(candidate_predictions["predictive_entropy_bits"].max()),
        "minimum_candidate_probability_margin": float(candidate_predictions["probability_margin"].min()),
        "contour_segment_count": int(contour["segment_id"].nunique()) if not contour.empty else 0,
        "probability_status": config["probability_status"],
        "fixed_test_labels_used": False,
        "exact_edge_orbit_used": bool(config["exact_edge_orbit_used"]),
        "held_out_CA_region_used": bool(config["held_out_CA_region_used"]),
        "random_seed_used": bool(config["random_seed_used"]),
        "claim_scope": config["claim_scope"],
    }

    dense_frame = pd.DataFrame(
        {
            "normalized_alpha_1": dense_X[:, 0],
            "normalized_alpha_2": dense_X[:, 1],
            "alpha_1": dense_X[:, 0] * alpha_scales[0],
            "alpha_2": dense_X[:, 1] * alpha_scales[1],
            "probability_C": dense_probability,
            "probability_B": 1.0 - dense_probability,
            "predictive_entropy_bits": predictive_entropy(dense_probability),
            "latent_mean": dense_mean,
            "latent_variance": dense_variance,
        }
    )

    outdir.mkdir(parents=True)
    shutil.copy2(CONFIG_FILE, outdir / CONFIG_FILE.name)
    pd.DataFrame(optimization_records).to_csv(outdir / "hyperparameter_optimization.csv", index=False)
    training_predictions.to_csv(outdir / "training_posterior_predictions.csv", index=False)
    loo.to_csv(outdir / "nested_LOO_predictions.csv", index=False)
    candidate_predictions.to_csv(outdir / "acquisition_pool_predictions.csv", index=False)
    active.to_csv(outdir / "active_acquisition_proposals.csv", index=False)
    baseline.to_csv(outdir / "space_filling_baseline_proposals.csv", index=False)
    fixed_test_predictions.to_csv(outdir / "fixed_test_blind_predictions.csv", index=False)
    dense_frame.to_csv(outdir / "dense_probability_grid.csv", index=False)
    contour.to_csv(outdir / "probability_half_contour.csv", index=False)
    np.savez_compressed(
        outdir / "probit_gp_model.npz",
        X=state.X,
        y=state.y,
        lengthscales=state.lengthscales,
        signal_std=np.array([state.signal_std]),
        jitter=np.array([state.jitter]),
        f=state.f,
        a=state.a,
        sqrt_W=state.sqrt_W,
        L=state.L,
    )
    with (outdir / "model_summary.json").open("w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2)
    metadata = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "runtime_seconds": time.perf_counter() - started,
        "python_version": sys.version,
        "numpy_version": np.__version__,
        "pandas_version": pd.__version__,
        "matplotlib_version": matplotlib.__version__,
        "working_directory": str(Path.cwd()),
        "git_commit": git_commit(),
    }
    with (outdir / "run_metadata.json").open("w", encoding="utf-8") as handle:
        json.dump(metadata, handle, indent=2)
    make_figure(
        outdir / "probabilistic_boundary_and_acquisition.png",
        training_predictions,
        fixed_test_predictions,
        active,
        baseline,
        grid_1,
        grid_2,
        probability_grid,
        entropy_grid,
        alpha_scales,
    )
    print(json.dumps(summary, indent=2))
    print(f"Saved probabilistic boundary model in: {outdir}")


if __name__ == "__main__":
    main()
