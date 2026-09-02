#!/usr/bin/env python3
"""Fit 12 matched sparse GPs and freeze four pairs for each method.

Only the 25 labels belonging to the current replicate are exposed during that
replicate. Candidate coordinates are taken from the geometry manifest, which
contains no outcomes. Evaluation labels are intentionally opened by a separate
script after this proposal directory has been hash-frozen.
"""

from __future__ import annotations

import hashlib
import json
import math
import time
from datetime import datetime, timezone
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist

from probit_gp import optimize_hyperparameters, predict_latent_joint, predict_proba


ROOT = Path(__file__).resolve().parent
CONFIG_FILE = ROOT / "ml_bc_matched_benchmark_protocol_config.json"
PROTOCOL_DIR = ROOT / "experiment_outputs" / "ml_BC_matched_benchmark_protocol_v1"
BANK_DIR = ROOT / "experiment_outputs" / "ml_BC_matched_label_bank_v1"
OUTDIR = ROOT / "experiment_outputs" / "ml_BC_matched_acquisition_proposals_v1"
METHODS = ("space_filling", "pointwise_uncertainty", "pair_aware")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def entropy_bits(probability: np.ndarray) -> np.ndarray:
    probability = np.clip(probability, 1e-12, 1.0 - 1e-12)
    return -(probability * np.log2(probability) + (1.0 - probability) * np.log2(1.0 - probability))


def candidate_pairs(candidate: pd.DataFrame, config: dict) -> pd.DataFrame:
    features = config["model"]["feature_columns"]
    X = candidate[features].to_numpy(float)
    distance = cdist(X, X)
    pair_cfg = config["pair_geometry"]
    lower = float(pair_cfg["minimum_normalized_separation"])
    upper = float(pair_cfg["maximum_normalized_separation"])
    first, second = np.where(np.triu((distance >= lower) & (distance <= upper), k=1))
    pairs = pd.DataFrame(
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
    low = np.minimum(pairs.first_state_index, pairs.second_state_index)
    high = np.maximum(pairs.first_state_index, pairs.second_state_index)
    pairs.insert(0, "pair_id", [f"pair_{a:03d}_{b:03d}" for a, b in zip(low, high)])
    return pairs


def common_pair_geometry(state, candidate: pd.DataFrame, pairs: pd.DataFrame, X_train, config: dict):
    features = config["model"]["feature_columns"]
    X = candidate[features].to_numpy(float)
    first = pairs.first_candidate_position.to_numpy(int)
    second = pairs.second_candidate_position.to_numpy(int)
    probability = predict_proba(state, X)[1]
    entropy = entropy_bits(probability)
    midpoint = 0.5 * (X[first] + X[second])
    direction = (X[second] - X[first]) / pairs.normalized_separation.to_numpy(float)[:, None]
    step = float(config["pair_geometry"]["normal_gradient_step"])
    plus_1 = predict_proba(state, midpoint + np.array([step, 0.0]))[1]
    minus_1 = predict_proba(state, midpoint - np.array([step, 0.0]))[1]
    plus_2 = predict_proba(state, midpoint + np.array([0.0, step]))[1]
    minus_2 = predict_proba(state, midpoint - np.array([0.0, step]))[1]
    gradient = np.column_stack([(plus_1 - minus_1) / (2 * step), (plus_2 - minus_2) / (2 * step)])
    gradient_norm = np.linalg.norm(gradient, axis=1)
    unit_normal = gradient / np.maximum(gradient_norm[:, None], 1e-14)
    raw_orientation = np.abs(np.sum(unit_normal * direction, axis=1))
    floor = float(config["pair_geometry"]["orientation_floor"])
    orientation = floor + (1.0 - floor) * raw_orientation
    scale = float(config["pair_geometry"]["distance_scale"])
    distance_penalty = np.exp(-0.5 * (pairs.normalized_separation.to_numpy(float) / scale) ** 2)
    novelty = cdist(X, X_train).min(axis=1)
    return {
        "X": X,
        "first": first,
        "second": second,
        "probability": probability,
        "entropy": entropy,
        "midpoint": midpoint,
        "gradient_norm": gradient_norm,
        "raw_orientation": raw_orientation,
        "orientation": orientation,
        "distance_penalty": distance_penalty,
        "novelty": novelty,
    }


def joint_opposite_probability(state, X: np.ndarray, first, second, sample_count, seed, eigenvalue_floor):
    mean, covariance = predict_latent_joint(state, X)
    eigenvalues, eigenvectors = np.linalg.eigh(covariance)
    clipped = np.maximum(eigenvalues, float(eigenvalue_floor))
    rng = np.random.default_rng(seed)
    standard = rng.standard_normal((int(sample_count), len(X)))
    # Apple's BLAS can emit spurious overflow warnings for this finite product;
    # the explicit finiteness gate below remains the operative numerical check.
    with np.errstate(over="ignore", divide="ignore", invalid="ignore"):
        samples = mean + standard @ (eigenvectors * np.sqrt(clipped)).T
    if not np.all(np.isfinite(samples)):
        raise FloatingPointError("Nonfinite joint posterior samples")
    opposite = np.mean((samples[:, first] >= 0.0) != (samples[:, second] >= 0.0), axis=0)
    return opposite, int(np.sum(eigenvalues < 0.0)), float(np.min(eigenvalues))


def rank_pairs(method: str, state, candidate, pairs, X_train, config, replicate):
    common = common_pair_geometry(state, candidate, pairs, X_train, config)
    first, second = common["first"], common["second"]
    ranked = pairs.copy()
    ranked["first_probability_C"] = common["probability"][first]
    ranked["second_probability_C"] = common["probability"][second]
    ranked["first_entropy_bits"] = common["entropy"][first]
    ranked["second_entropy_bits"] = common["entropy"][second]
    ranked["first_novelty"] = common["novelty"][first]
    ranked["second_novelty"] = common["novelty"][second]
    ranked["distance_penalty"] = common["distance_penalty"]
    ranked["raw_normal_orientation"] = common["raw_orientation"]
    ranked["orientation_factor"] = common["orientation"]
    ranked["midpoint_gradient_norm"] = common["gradient_norm"]
    diagnostics = {}
    if method == "space_filling":
        information = np.sqrt(common["novelty"][first] * common["novelty"][second])
        score = information * common["distance_penalty"]
        ranked["information_factor"] = information
        ranked["joint_opposite_sign_probability"] = np.nan
    elif method == "pointwise_uncertainty":
        information = np.sqrt(common["entropy"][first] * common["entropy"][second])
        score = information * common["distance_penalty"] * common["orientation"]
        ranked["information_factor"] = information
        ranked["joint_opposite_sign_probability"] = np.nan
    elif method == "pair_aware":
        sampling = config["joint_posterior_sampling"]
        opposite, negative_count, minimum_eigenvalue = joint_opposite_probability(
            state,
            common["X"],
            first,
            second,
            int(sampling["sample_count"]),
            int(sampling["base_random_seed"]) + replicate - 1,
            float(sampling["eigenvalue_floor"]),
        )
        score = opposite * common["distance_penalty"] * common["orientation"]
        ranked["information_factor"] = opposite
        ranked["joint_opposite_sign_probability"] = opposite
        diagnostics = {
            "negative_joint_covariance_eigenvalue_count_before_clipping": negative_count,
            "minimum_joint_covariance_eigenvalue_before_clipping": minimum_eigenvalue,
        }
    else:
        raise ValueError(f"Unknown acquisition method {method}")
    ranked["acquisition_score"] = score
    ranked = ranked.sort_values(
        ["acquisition_score", "information_factor", "normalized_separation", "pair_id"],
        ascending=[False, False, True, True],
    ).reset_index(drop=True)
    ranked.insert(1, "score_rank", np.arange(1, len(ranked) + 1))
    return ranked, diagnostics


def greedy_batch(ranking: pd.DataFrame, config: dict) -> pd.DataFrame:
    batch_cfg = config["batch"]
    pair_count = int(batch_cfg["pair_count"])
    diversity_scale = float(batch_cfg["midpoint_diversity_scale"])
    diversity_floor = float(batch_cfg["diversity_floor"])
    available = ranking.copy()
    available["midpoint_1"] = 0.5 * (
        available.first_normalized_alpha_1 + available.second_normalized_alpha_1
    )
    available["midpoint_2"] = 0.5 * (
        available.first_normalized_alpha_2 + available.second_normalized_alpha_2
    )
    selected = []
    used_endpoints: set[int] = set()
    for step in range(1, pair_count + 1):
        choices = available[
            ~available.first_state_index.isin(used_endpoints)
            & ~available.second_state_index.isin(used_endpoints)
        ].copy()
        if choices.empty:
            raise ValueError("No disjoint pair remains for the declared eight-state budget")
        if not selected:
            choices["batch_diversity_factor"] = 1.0
        else:
            midpoint = choices[["midpoint_1", "midpoint_2"]].to_numpy(float)
            chosen_midpoint = np.array([[row["midpoint_1"], row["midpoint_2"]] for row in selected])
            distance = cdist(midpoint, chosen_midpoint).min(axis=1)
            choices["batch_diversity_factor"] = diversity_floor + (1.0 - diversity_floor) * np.minimum(
                distance / diversity_scale, 1.0
            )
        choices["batch_score"] = choices.acquisition_score * choices.batch_diversity_factor
        choices = choices.sort_values(
            ["batch_score", "acquisition_score", "pair_id"], ascending=[False, False, True]
        )
        chosen = choices.iloc[0].to_dict()
        chosen["proposal_rank"] = step
        selected.append(chosen)
        used_endpoints |= {int(chosen["first_state_index"]), int(chosen["second_state_index"])}
    columns = ["proposal_rank"] + [item for item in selected[0] if item != "proposal_rank"]
    return pd.DataFrame(selected)[columns]


def make_figure(path: Path, geometry, designs, proposals):
    figure, axes = plt.subplots(3, 4, figsize=(12, 9), constrained_layout=True)
    colors = {"space_filling": "#009E73", "pointwise_uncertainty": "#E69F00", "pair_aware": "#CC79A7"}
    for replicate, axis in enumerate(axes.flat, start=1):
        axis.scatter(geometry.normalized_alpha_1, geometry.normalized_alpha_2, s=6, color="0.88")
        design = designs[designs.replicate == replicate]
        axis.scatter(design.normalized_alpha_1, design.normalized_alpha_2, s=13, color="#0072B2")
        for method in METHODS:
            batch = proposals[(proposals.replicate == replicate) & (proposals.method == method)]
            for row in batch.itertuples(index=False):
                axis.plot(
                    [row.first_normalized_alpha_1, row.second_normalized_alpha_1],
                    [row.first_normalized_alpha_2, row.second_normalized_alpha_2],
                    "-o",
                    color=colors[method],
                    linewidth=1.0,
                    markersize=2.5,
                    alpha=0.85,
                )
        axis.set_title(f"design {replicate:02d}")
        axis.set_aspect("equal")
        axis.set_xlim(-1.05, 1.05)
        axis.set_ylim(-1.05, 1.05)
    figure.suptitle("Frozen matched four-pair batches; candidate outcomes unopened")
    figure.supxlabel(r"normalized $\alpha_1$")
    figure.supylabel(r"normalized $\alpha_2$")
    figure.savefig(path, dpi=250)
    plt.close(figure)


def main() -> None:
    started = time.perf_counter()
    if OUTDIR.exists():
        raise FileExistsError(f"Refusing to overwrite frozen proposals: {OUTDIR}")
    config = json.loads(CONFIG_FILE.read_text(encoding="utf-8"))
    protocol = json.loads((PROTOCOL_DIR / "frozen_protocol_manifest.json").read_text(encoding="utf-8"))
    bank_manifest_path = BANK_DIR / "complete_label_bank_manifest.json"
    bank_manifest = json.loads(bank_manifest_path.read_text(encoding="utf-8"))
    if not protocol.get("protocol_frozen") or not bank_manifest.get("binary_BC_plane_gate_passed"):
        raise ValueError("Frozen protocol or controlled binary-plane gate failed")
    bank_path = BANK_DIR / "complete_label_bank.csv"
    if sha256(bank_path) != bank_manifest["complete_label_bank_sha256"]:
        raise ValueError("Complete label bank changed")

    geometry = pd.read_csv(ROOT / config["geometry_manifest_csv"]).sort_values("state_index")
    bank = pd.read_csv(bank_path, usecols=["state_index", "final_label"])
    designs_path = PROTOCOL_DIR / "matched_initial_designs.csv"
    if sha256(designs_path) != protocol["matched_initial_designs_sha256"]:
        raise ValueError("Matched initial designs changed")
    designs = pd.read_csv(designs_path)
    model_cfg = config["model"]
    features = model_cfg["feature_columns"]

    proposal_frames = []
    training_frames = []
    hyperparameter_rows = []
    diagnostic_rows = []
    for replicate in range(1, int(config["initial_design"]["replicate_count"]) + 1):
        design = designs[designs.replicate == replicate].copy()
        training = design.merge(bank, on="state_index", how="left", validate="one_to_one")
        if len(training) != int(config["initial_design"]["training_count"]):
            raise ValueError(f"Replicate {replicate} training count changed")
        if set(training.final_label) != {"B", "C"}:
            raise ValueError(f"Replicate {replicate} does not contain both B and C only")
        training["label_exposed_during_fit"] = True
        training_frames.append(training)
        X_train = training[features].to_numpy(float)
        y_train = np.where(training.final_label == model_cfg["positive_class"], 1.0, -1.0)
        state, optimization = optimize_hyperparameters(
            X_train,
            y_train,
            model_cfg["initial_starts"],
            lengthscale_bounds=tuple(model_cfg["lengthscale_bounds"]),
            signal_std_bounds=tuple(model_cfg["signal_std_bounds"]),
            jitter=float(model_cfg["jitter"]),
            maximum_laplace_iterations=int(model_cfg["maximum_laplace_iterations"]),
            laplace_tolerance=float(model_cfg["laplace_tolerance"]),
            maximum_optimizer_iterations=int(model_cfg["maximum_optimizer_iterations"]),
            function_tolerance=float(model_cfg["function_tolerance"]),
            log_prior_mean=np.log(np.asarray(model_cfg["log_prior_mean_parameters"], dtype=float)),
            log_prior_std=np.asarray(model_cfg["log_prior_std"], dtype=float),
        )
        if not state.converged:
            raise RuntimeError(f"GP Laplace fit failed for replicate {replicate}")
        hyperparameter_rows.append(
            {
                "replicate": replicate,
                "B_training_count": int((training.final_label == "B").sum()),
                "C_training_count": int((training.final_label == "C").sum()),
                "lengthscale_1": state.lengthscales[0],
                "lengthscale_2": state.lengthscales[1],
                "signal_std": state.signal_std,
                "log_posterior_objective": state.log_marginal_likelihood,
                "optimizer_start_count": len(model_cfg["initial_starts"]),
                "optimizer_record_count": len(optimization),
            }
        )
        candidate = geometry[~geometry.state_index.isin(training.state_index)].copy().reset_index(drop=True)
        if "final_label" in candidate.columns:
            raise RuntimeError("Candidate outcomes leaked into the acquisition frame")
        pairs = candidate_pairs(candidate, config)
        for method in METHODS:
            ranking, diagnostics = rank_pairs(method, state, candidate, pairs, X_train, config, replicate)
            batch = greedy_batch(ranking, config)
            batch.insert(0, "method", method)
            batch.insert(0, "replicate", replicate)
            proposal_frames.append(batch)
            diagnostic_rows.append({"replicate": replicate, "method": method, **diagnostics})
        print(f"replicate {replicate:02d}: frozen {len(METHODS) * int(config['batch']['pair_count'])} pairs")

    proposals = pd.concat(proposal_frames, ignore_index=True)
    training_labels = pd.concat(training_frames, ignore_index=True)
    expected = int(config["initial_design"]["replicate_count"]) * len(METHODS) * int(config["batch"]["pair_count"])
    if len(proposals) != expected:
        raise ValueError("Frozen proposal count changed")
    endpoint_count = proposals.groupby(["replicate", "method"]).apply(
        lambda frame: len(set(frame.first_state_index) | set(frame.second_state_index))
    )
    if not np.all(endpoint_count == int(config["batch"]["new_PDE_endpoint_budget"])):
        raise ValueError("A method does not use exactly eight disjoint endpoints")

    OUTDIR.mkdir(parents=True)
    proposal_path = OUTDIR / "frozen_matched_pair_proposals.csv"
    training_path = OUTDIR / "exposed_training_labels.csv"
    hyperparameter_path = OUTDIR / "replicate_GP_hyperparameters.csv"
    diagnostics_path = OUTDIR / "acquisition_numerical_diagnostics.csv"
    proposals.to_csv(proposal_path, index=False)
    training_labels.to_csv(training_path, index=False)
    pd.DataFrame(hyperparameter_rows).to_csv(hyperparameter_path, index=False)
    pd.DataFrame(diagnostic_rows).to_csv(diagnostics_path, index=False)
    make_figure(OUTDIR / "frozen_matched_pair_proposals.png", geometry, designs, proposals)
    manifest = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "proposal_batches_frozen": True,
        "replicate_count": int(config["initial_design"]["replicate_count"]),
        "methods": list(METHODS),
        "pair_count_per_method_per_replicate": int(config["batch"]["pair_count"]),
        "endpoint_budget_per_method_per_replicate": int(config["batch"]["new_PDE_endpoint_budget"]),
        "proposal_count": len(proposals),
        "candidate_outcomes_used": False,
        "only_current_replicate_training_labels_exposed": True,
        "exact_edge_orbit_used": False,
        "heldout_CBA_region_used": False,
        "frozen_matched_pair_proposals_sha256": sha256(proposal_path),
        "exposed_training_labels_sha256": sha256(training_path),
        "replicate_GP_hyperparameters_sha256": sha256(hyperparameter_path),
        "complete_label_bank_sha256": bank_manifest["complete_label_bank_sha256"],
        "frozen_protocol_manifest_sha256": sha256(PROTOCOL_DIR / "frozen_protocol_manifest.json"),
        "acquisition_script_sha256": sha256(Path(__file__)),
        "runtime_seconds": time.perf_counter() - started,
        "interpretation": config["interpretation"],
    }
    (OUTDIR / "frozen_proposal_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(json.dumps(manifest, indent=2))
    print(f"Frozen matched proposals saved in: {OUTDIR}")


if __name__ == "__main__":
    main()
