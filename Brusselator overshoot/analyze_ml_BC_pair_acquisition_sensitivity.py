#!/usr/bin/env python3
"""Develop and stress-test joint-posterior B--C pair acquisition.

Only the original 25 labels enter model fitting. Previously opened outcomes
are read solely as state indices so that new proposals cannot reuse them.
This stage does not designate a primary batch and does not run the PDE.
"""

import json
import math
import os
import shutil
import time
from datetime import datetime, timezone
from pathlib import Path
from types import SimpleNamespace

os.environ.setdefault("MPLCONFIGDIR", "/private/tmp/brusselator_matplotlib")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist
from scipy.stats import spearmanr

from probit_gp import (
    LaplaceState,
    ard_rbf_kernel,
    optimize_hyperparameters,
    predict_latent_joint,
    predict_proba,
)


ROOT = Path(__file__).resolve().parent
CONFIG_FILE = ROOT / "ml_bc_pair_acquisition_sensitivity_config.json"


def load_saved_state(path):
    data = np.load(path)
    X = data["X"]
    lengthscales = data["lengthscales"]
    signal_std = float(data["signal_std"][0])
    jitter = float(data["jitter"][0])
    K = ard_rbf_kernel(X, X, lengthscales, signal_std) + jitter * np.eye(len(X))
    return LaplaceState(
        X=X,
        y=data["y"],
        lengthscales=lengthscales,
        signal_std=signal_std,
        jitter=jitter,
        K=K,
        f=data["f"],
        a=data["a"],
        sqrt_W=data["sqrt_W"],
        L=data["L"],
        log_marginal_likelihood=np.nan,
        iterations=0,
        converged=True,
    )


def fit_sensitivity_model(X, y, fit, config):
    if fit["kind"] == "saved_pretest_model":
        return load_saved_state(ROOT / config["nominal_model_npz"]), pd.DataFrame()
    optimization = config["optimization"]
    kwargs = {}
    if fit["kind"] == "refit_with_log_parameter_prior":
        kwargs["log_prior_mean"] = np.log(np.asarray(fit["log_prior_mean_parameters"], dtype=float))
        kwargs["log_prior_std"] = np.asarray(fit["log_prior_std"], dtype=float)
    state, records = optimize_hyperparameters(
        X,
        y,
        optimization["initial_starts"],
        lengthscale_bounds=tuple(fit["lengthscale_bounds"]),
        signal_std_bounds=tuple(fit["signal_std_bounds"]),
        jitter=float(optimization["jitter"]),
        maximum_laplace_iterations=int(optimization["maximum_laplace_iterations"]),
        laplace_tolerance=float(optimization["laplace_tolerance"]),
        maximum_optimizer_iterations=int(optimization["maximum_optimizer_iterations"]),
        function_tolerance=float(optimization["function_tolerance"]),
        **kwargs,
    )
    if not state.converged:
        raise RuntimeError(f"Laplace iteration failed for {fit['name']}")
    return state, pd.DataFrame(records)


def candidate_pairs(candidate, geometry):
    features = ["normalized_alpha_1", "normalized_alpha_2"]
    X = candidate[features].to_numpy(float)
    distance = cdist(X, X)
    lower = float(geometry["minimum_normalized_separation"])
    upper = float(geometry["maximum_normalized_separation"])
    first, second = np.where(np.triu((distance >= lower) & (distance <= upper), k=1))
    rows = pd.DataFrame(
        {
            "first_candidate_position": first,
            "second_candidate_position": second,
            "first_state_index": candidate.iloc[first].state_index.to_numpy(int),
            "second_state_index": candidate.iloc[second].state_index.to_numpy(int),
            "first_sample_id": candidate.iloc[first].sample_id.to_numpy(),
            "second_sample_id": candidate.iloc[second].sample_id.to_numpy(),
            "first_normalized_alpha_1": X[first, 0],
            "first_normalized_alpha_2": X[first, 1],
            "second_normalized_alpha_1": X[second, 0],
            "second_normalized_alpha_2": X[second, 1],
            "normalized_separation": distance[first, second],
        }
    )
    low = np.minimum(rows.first_state_index, rows.second_state_index)
    high = np.maximum(rows.first_state_index, rows.second_state_index)
    rows.insert(0, "pair_id", [f"pair_{a:03d}_{b:03d}" for a, b in zip(low, high)])
    return rows


def joint_latent_samples(state, X, sampling):
    mean, covariance = predict_latent_joint(state, X)
    eigenvalues, eigenvectors = np.linalg.eigh(covariance)
    clipped = np.maximum(eigenvalues, float(sampling["eigenvalue_floor"]))
    rng = np.random.default_rng(int(sampling["random_seed"]))
    standard = rng.standard_normal((int(sampling["sample_count"]), len(X)))
    samples = mean + standard @ (eigenvectors * np.sqrt(clipped)).T
    return mean, covariance, samples, int(np.sum(eigenvalues < 0.0)), float(np.min(eigenvalues))


def score_pairs(state, candidate, pairs, config):
    features = config["feature_columns"]
    X = candidate[features].to_numpy(float)
    mean, covariance, samples, negative_count, minimum_eigenvalue = joint_latent_samples(
        state, X, config["joint_posterior_sampling"]
    )
    first = pairs.first_candidate_position.to_numpy(int)
    second = pairs.second_candidate_position.to_numpy(int)
    opposite = np.mean((samples[:, first] >= 0.0) != (samples[:, second] >= 0.0), axis=0)

    first_X = X[first]
    second_X = X[second]
    midpoint = 0.5 * (first_X + second_X)
    direction = (second_X - first_X) / pairs.normalized_separation.to_numpy(float)[:, None]
    step = float(config["pair_geometry"]["normal_gradient_step"])
    offsets_1 = np.array([step, 0.0])
    offsets_2 = np.array([0.0, step])
    plus_1 = predict_proba(state, midpoint + offsets_1)[1]
    minus_1 = predict_proba(state, midpoint - offsets_1)[1]
    plus_2 = predict_proba(state, midpoint + offsets_2)[1]
    minus_2 = predict_proba(state, midpoint - offsets_2)[1]
    gradient = np.column_stack([(plus_1 - minus_1) / (2 * step), (plus_2 - minus_2) / (2 * step)])
    gradient_norm = np.linalg.norm(gradient, axis=1)
    unit_normal = gradient / np.maximum(gradient_norm[:, None], 1e-14)
    raw_orientation = np.abs(np.sum(unit_normal * direction, axis=1))
    floor = float(config["pair_geometry"]["orientation_floor"])
    orientation = floor + (1.0 - floor) * raw_orientation
    scale = float(config["pair_geometry"]["distance_scale"])
    distance_penalty = np.exp(-0.5 * (pairs.normalized_separation.to_numpy(float) / scale) ** 2)
    score = opposite * distance_penalty * orientation

    result = pairs.copy()
    result["joint_opposite_sign_probability"] = opposite
    result["distance_penalty"] = distance_penalty
    result["raw_normal_orientation"] = raw_orientation
    result["orientation_factor"] = orientation
    result["midpoint_gradient_norm"] = gradient_norm
    result["acquisition_score"] = score
    result["first_latent_mean"] = mean[first]
    result["second_latent_mean"] = mean[second]
    result["first_latent_variance"] = np.diag(covariance)[first]
    result["second_latent_variance"] = np.diag(covariance)[second]
    result = result.sort_values(
        ["acquisition_score", "joint_opposite_sign_probability", "normalized_separation", "pair_id"],
        ascending=[False, False, True, True],
    ).reset_index(drop=True)
    result.insert(1, "score_rank", np.arange(1, len(result) + 1))
    diagnostics = {
        "negative_joint_covariance_eigenvalue_count_before_clipping": negative_count,
        "minimum_joint_covariance_eigenvalue_before_clipping": minimum_eigenvalue,
        "maximum_score": float(result.acquisition_score.max()),
        "maximum_joint_opposite_sign_probability": float(result.joint_opposite_sign_probability.max()),
    }
    return result, diagnostics


def greedy_batch(ranking, config):
    pair_count = int(config["batch"]["pair_count"])
    diversity_scale = float(config["batch"]["midpoint_diversity_scale"])
    diversity_floor = float(config["batch"]["diversity_floor"])
    available = ranking.copy()
    available["midpoint_1"] = 0.5 * (
        available.first_normalized_alpha_1 + available.second_normalized_alpha_1
    )
    available["midpoint_2"] = 0.5 * (
        available.first_normalized_alpha_2 + available.second_normalized_alpha_2
    )
    selected = []
    used_endpoints = set()
    for step in range(1, pair_count + 1):
        candidates = available[
            ~available.first_state_index.isin(used_endpoints)
            & ~available.second_state_index.isin(used_endpoints)
        ].copy()
        if not selected:
            candidates["batch_diversity_factor"] = 1.0
        else:
            midpoints = candidates[["midpoint_1", "midpoint_2"]].to_numpy(float)
            selected_midpoints = np.asarray(
                [[row["midpoint_1"], row["midpoint_2"]] for row in selected]
            )
            distance = np.min(cdist(midpoints, selected_midpoints), axis=1)
            candidates["batch_diversity_factor"] = diversity_floor + (1.0 - diversity_floor) * np.minimum(
                distance / diversity_scale, 1.0
            )
        candidates["batch_score"] = candidates.acquisition_score * candidates.batch_diversity_factor
        candidates = candidates.sort_values(
            ["batch_score", "acquisition_score", "pair_id"], ascending=[False, False, True]
        )
        chosen = candidates.iloc[0].to_dict()
        chosen["batch_step"] = step
        selected.append(chosen)
        used_endpoints |= {int(chosen["first_state_index"]), int(chosen["second_state_index"])}
    columns = ["batch_step"] + [column for column in selected[0] if column != "batch_step"]
    return pd.DataFrame(selected)[columns]


def ranking_comparisons(rankings, batches, top_k):
    names = list(rankings)
    rows = []
    for left_index, left in enumerate(names):
        for right in names[left_index + 1 :]:
            left_rank = rankings[left].set_index("pair_id")
            right_rank = rankings[right].set_index("pair_id")
            common = left_rank.index.intersection(right_rank.index)
            correlation = spearmanr(
                left_rank.loc[common, "score_rank"], right_rank.loc[common, "score_rank"]
            ).statistic
            left_top = set(left_rank.nsmallest(top_k, "score_rank").index)
            right_top = set(right_rank.nsmallest(top_k, "score_rank").index)
            left_batch = set(batches[left].pair_id)
            right_batch = set(batches[right].pair_id)
            rows.append(
                {
                    "left_fit": left,
                    "right_fit": right,
                    "all_pair_spearman_rank_correlation": float(correlation),
                    "top_k": top_k,
                    "top_k_intersection": len(left_top & right_top),
                    "top_k_jaccard": len(left_top & right_top) / len(left_top | right_top),
                    "four_pair_batch_intersection": len(left_batch & right_batch),
                    "four_pair_batch_jaccard": len(left_batch & right_batch) / len(left_batch | right_batch),
                }
            )
    return pd.DataFrame(rows)


def make_figure(path, candidate, rankings, batches):
    names = list(rankings)
    figure, axes = plt.subplots(1, len(names), figsize=(5.1 * len(names), 4.7), constrained_layout=True)
    for axis, name in zip(np.atleast_1d(axes), names):
        ranking = rankings[name]
        image = axis.scatter(
            0.5 * (ranking.first_normalized_alpha_1 + ranking.second_normalized_alpha_1),
            0.5 * (ranking.first_normalized_alpha_2 + ranking.second_normalized_alpha_2),
            c=ranking.acquisition_score,
            s=12,
            cmap="viridis",
            alpha=0.7,
        )
        for row in batches[name].itertuples(index=False):
            axis.plot(
                [row.first_normalized_alpha_1, row.second_normalized_alpha_1],
                [row.first_normalized_alpha_2, row.second_normalized_alpha_2],
                "-o",
                color="magenta",
                linewidth=2,
                markersize=4,
            )
            axis.text(row.midpoint_1, row.midpoint_2, str(int(row.batch_step)), fontsize=8)
        axis.scatter(
            candidate.normalized_alpha_1,
            candidate.normalized_alpha_2,
            facecolors="none",
            edgecolors="0.6",
            s=12,
            linewidths=0.4,
        )
        axis.set_title(name)
        axis.set_xlabel(r"normalized $\alpha_1$")
        axis.set_ylabel(r"normalized $\alpha_2$")
        axis.set_aspect("equal")
        figure.colorbar(image, ax=axis, label="pair acquisition score")
    figure.suptitle("Joint-posterior pair acquisition sensitivity; no endpoint outcomes used")
    figure.savefig(path, dpi=300)
    plt.close(figure)


def main():
    started = time.perf_counter()
    with CONFIG_FILE.open("r", encoding="utf-8") as handle:
        config = json.load(handle)
    outdir = ROOT / config["output_root"] / config["experiment_name"]
    if outdir.exists():
        raise FileExistsError(f"Refusing to overwrite sensitivity output: {outdir}")

    manifest = pd.read_csv(ROOT / config["geometry_manifest_csv"])
    labels = pd.read_csv(ROOT / config["initial_labels_csv"])
    training = manifest.merge(
        labels[["state_index", "final_label"]], on="state_index", validate="one_to_one"
    ).sort_values("state_index")
    if len(training) != config["training_count"] or not np.all(training.design_role == config["training_role"]):
        raise ValueError("The original 25-state training set changed")
    if set(training.final_label) - {"B", "C"}:
        raise ValueError("Pair acquisition requires binary B/C training labels")

    revealed = set()
    revealed_counts = {}
    for relative in config["revealed_index_sources"]:
        frame = pd.read_csv(ROOT / relative, usecols=["state_index"])
        values = set(frame.state_index.astype(int))
        revealed |= values
        revealed_counts[relative] = len(values)
    candidate = manifest[
        (manifest.design_role == config["candidate_role"])
        & ~manifest.state_index.isin(revealed)
    ].copy().sort_values("state_index").reset_index(drop=True)
    if candidate.empty:
        raise ValueError("No unopened acquisition-pool states remain")
    pairs = candidate_pairs(candidate, config["pair_geometry"])
    if pairs.empty:
        raise ValueError("No admissible nearby pairs remain")

    features = config["feature_columns"]
    X_train = training[features].to_numpy(float)
    y_train = np.where(training.final_label == config["positive_class"], 1.0, -1.0)
    rankings = {}
    batches = {}
    diagnostics = {}
    hyperparameter_rows = []
    optimization_tables = {}
    for fit in config["ranking_sensitivity"]["fits"]:
        state, optimization = fit_sensitivity_model(X_train, y_train, fit, config)
        ranking, joint_diagnostics = score_pairs(state, candidate, pairs, config)
        batch = greedy_batch(ranking, config)
        rankings[fit["name"]] = ranking
        batches[fit["name"]] = batch
        diagnostics[fit["name"]] = joint_diagnostics
        optimization_tables[fit["name"]] = optimization
        bounds = fit.get("lengthscale_bounds", [np.nan, np.nan])
        signal_bounds = fit.get("signal_std_bounds", [np.nan, np.nan])
        hyperparameter_rows.append(
            {
                "fit": fit["name"],
                "kind": fit["kind"],
                "lengthscale_1": state.lengthscales[0],
                "lengthscale_2": state.lengthscales[1],
                "signal_std": state.signal_std,
                "lengthscale_at_upper_bound": bool(
                    np.isfinite(bounds[1]) and np.any(state.lengthscales >= 0.999999 * bounds[1])
                ),
                "signal_std_at_upper_bound": bool(
                    np.isfinite(signal_bounds[1]) and state.signal_std >= 0.999999 * signal_bounds[1]
                ),
                "log_marginal_likelihood": state.log_marginal_likelihood,
            }
        )

    comparisons = ranking_comparisons(
        rankings, batches, int(config["ranking_sensitivity"]["top_k"])
    )
    hyperparameters = pd.DataFrame(hyperparameter_rows)
    summary = {
        "experiment_name": config["experiment_name"],
        "phase": config["phase"],
        "training_count": len(training),
        "candidate_endpoint_count": len(candidate),
        "admissible_pair_count": len(pairs),
        "revealed_state_count_excluded": len(revealed),
        "revealed_index_source_counts": revealed_counts,
        "fit_count": len(rankings),
        "minimum_top_k_jaccard": float(comparisons.top_k_jaccard.min()),
        "minimum_all_pair_spearman_rank_correlation": float(
            comparisons.all_pair_spearman_rank_correlation.min()
        ),
        "minimum_four_pair_batch_intersection": int(comparisons.four_pair_batch_intersection.min()),
        "no_outcome_column_read_from_revealed_sources": True,
        "new_PDE_outcomes_used": False,
        "exact_edge_orbit_used": False,
        "held_out_CA_region_used": False,
        "primary_batch_frozen": False,
        "interpretation": "training-only method-development sensitivity; review ranking robustness before freezing an eight-run pair batch",
    }

    outdir.mkdir(parents=True)
    shutil.copy2(CONFIG_FILE, outdir / CONFIG_FILE.name)
    candidate.to_csv(outdir / "unopened_candidate_endpoints.csv", index=False)
    pairs.to_csv(outdir / "admissible_pair_geometry.csv", index=False)
    hyperparameters.to_csv(outdir / "hyperparameter_sensitivity.csv", index=False)
    comparisons.to_csv(outdir / "pair_ranking_sensitivity.csv", index=False)
    for name in rankings:
        rankings[name].to_csv(outdir / f"{name}_pair_ranking.csv", index=False)
        batches[name].to_csv(outdir / f"{name}_four_pair_batch.csv", index=False)
        optimization_tables[name].to_csv(outdir / f"{name}_optimization.csv", index=False)
    with (outdir / "joint_posterior_diagnostics.json").open("w", encoding="utf-8") as handle:
        json.dump(diagnostics, handle, indent=2)
    with (outdir / "sensitivity_summary.json").open("w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2)
    with (outdir / "run_metadata.json").open("w", encoding="utf-8") as handle:
        json.dump(
            {
                "created_utc": datetime.now(timezone.utc).isoformat(),
                "runtime_seconds": time.perf_counter() - started,
                "numpy_version": np.__version__,
                "pandas_version": pd.__version__,
                "matplotlib_version": matplotlib.__version__,
            },
            handle,
            indent=2,
        )
    make_figure(outdir / "pair_acquisition_sensitivity.png", candidate, rankings, batches)
    print(hyperparameters.to_string(index=False))
    print(comparisons.to_string(index=False))
    print(json.dumps(summary, indent=2))
    print(f"Pair-acquisition sensitivity saved in: {outdir}")


if __name__ == "__main__":
    main()
