#!/usr/bin/env python3
"""Freeze masked-replay proposals for the B--C acquisition robustness suite.

The full label-bank hash is verified, but only labels in the current matched
training prefix are merged into a GP fit. Candidate frames contain geometry
only. Outcome opening is delegated to a separate evaluation script after this
proposal table and its manifest have been written and hashed.
"""

from __future__ import annotations

import hashlib
import json
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist

from probit_gp_kernels import optimize_hyperparameters, predict_latent_joint, predict_proba


ROOT = Path(__file__).resolve().parent
CONFIG_FILE = ROOT / "ml_bc_acquisition_robustness_config.json"


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
    if pairs.empty:
        raise ValueError("No admissible candidate pairs remain")
    return pairs


def common_pair_geometry(state, candidate: pd.DataFrame, pairs: pd.DataFrame, X_train, config: dict):
    features = config["model"]["feature_columns"]
    X = candidate[features].to_numpy(float)
    first = pairs.first_candidate_position.to_numpy(int)
    second = pairs.second_candidate_position.to_numpy(int)
    probability = predict_proba(state, X)[0]
    entropy = entropy_bits(probability)
    midpoint = 0.5 * (X[first] + X[second])
    separation = pairs.normalized_separation.to_numpy(float)
    direction = (X[second] - X[first]) / separation[:, None]
    step = float(config["pair_geometry"]["normal_gradient_step"])
    plus_1 = predict_proba(state, midpoint + np.array([step, 0.0]))[0]
    minus_1 = predict_proba(state, midpoint - np.array([step, 0.0]))[0]
    plus_2 = predict_proba(state, midpoint + np.array([0.0, step]))[0]
    minus_2 = predict_proba(state, midpoint - np.array([0.0, step]))[0]
    gradient = np.column_stack([(plus_1 - minus_1) / (2 * step), (plus_2 - minus_2) / (2 * step)])
    gradient_norm = np.linalg.norm(gradient, axis=1)
    unit_normal = gradient / np.maximum(gradient_norm[:, None], 1e-14)
    raw_orientation = np.abs(np.sum(unit_normal * direction, axis=1))
    floor = float(config["pair_geometry"]["orientation_floor"])
    orientation = floor + (1.0 - floor) * raw_orientation
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
        "novelty": novelty,
        "separation": separation,
    }


def joint_opposite_probability(state, X, first, second, sampling: dict, seed: int):
    mean, covariance = predict_latent_joint(state, X)
    eigenvalues, eigenvectors = np.linalg.eigh(covariance)
    clipped = np.maximum(eigenvalues, float(sampling["eigenvalue_floor"]))
    rng = np.random.default_rng(seed)
    standard = rng.standard_normal((int(sampling["sample_count"]), len(X)))
    with np.errstate(over="ignore", divide="ignore", invalid="ignore"):
        samples = mean + standard @ (eigenvectors * np.sqrt(clipped)).T
    if not np.all(np.isfinite(samples)):
        raise FloatingPointError("Nonfinite joint posterior samples")
    opposite = np.mean((samples[:, first] >= 0.0) != (samples[:, second] >= 0.0), axis=0)
    return opposite, int(np.sum(eigenvalues < 0.0)), float(np.min(eigenvalues))


def score_ranking(
    method: str,
    pair_score: str,
    distance_scale: float,
    pairs: pd.DataFrame,
    common: dict,
    opposite_probability: np.ndarray,
) -> pd.DataFrame:
    first = common["first"]
    second = common["second"]
    distance_penalty = np.exp(-0.5 * (common["separation"] / distance_scale) ** 2)
    ranked = pairs.copy()
    ranked["first_probability_C"] = common["probability"][first]
    ranked["second_probability_C"] = common["probability"][second]
    ranked["first_entropy_bits"] = common["entropy"][first]
    ranked["second_entropy_bits"] = common["entropy"][second]
    ranked["first_novelty"] = common["novelty"][first]
    ranked["second_novelty"] = common["novelty"][second]
    ranked["joint_opposite_sign_probability"] = opposite_probability
    ranked["distance_penalty"] = distance_penalty
    ranked["raw_normal_orientation"] = common["raw_orientation"]
    ranked["orientation_factor"] = common["orientation"]
    ranked["midpoint_gradient_norm"] = common["gradient_norm"]

    if method == "space_filling":
        information = np.sqrt(common["novelty"][first] * common["novelty"][second])
        score = information * distance_penalty
        score_definition = "sqrt(endpoint_novelty_product)*D"
    elif method == "pointwise_uncertainty":
        information = np.sqrt(common["entropy"][first] * common["entropy"][second])
        score = information * distance_penalty * common["orientation"]
        score_definition = "sqrt(endpoint_entropy_product)*D*O"
    elif method == "pair_aware":
        information = opposite_probability
        if pair_score == "q":
            score = information
        elif pair_score == "qD":
            score = information * distance_penalty
        elif pair_score == "qDO":
            score = information * distance_penalty * common["orientation"]
        else:
            raise ValueError(f"Unknown pair-aware score: {pair_score}")
        score_definition = pair_score
    else:
        raise ValueError(f"Unknown method: {method}")
    ranked["information_factor"] = information
    ranked["acquisition_score"] = score
    ranked["score_definition"] = score_definition
    ranked = ranked.sort_values(
        ["acquisition_score", "information_factor", "normalized_separation", "pair_id"],
        ascending=[False, False, True, True],
    ).reset_index(drop=True)
    ranked.insert(1, "score_rank", np.arange(1, len(ranked) + 1))
    return ranked


def greedy_batch(ranking: pd.DataFrame, config: dict) -> pd.DataFrame:
    batch_cfg = config["batch"]
    pair_count = int(batch_cfg["pair_count"])
    diversity_scale = float(batch_cfg["midpoint_diversity_scale"])
    diversity_floor = float(batch_cfg["diversity_floor"])
    available = ranking.copy()
    available["midpoint_1"] = 0.5 * (available.first_normalized_alpha_1 + available.second_normalized_alpha_1)
    available["midpoint_2"] = 0.5 * (available.first_normalized_alpha_2 + available.second_normalized_alpha_2)
    selected = []
    used_endpoints: set[int] = set()
    for step in range(1, pair_count + 1):
        choices = available[
            ~available.first_state_index.isin(used_endpoints)
            & ~available.second_state_index.isin(used_endpoints)
        ].copy()
        if choices.empty:
            raise ValueError("No disjoint pair remains for the declared endpoint budget")
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
    columns = ["proposal_rank"] + [column for column in selected[0] if column != "proposal_rank"]
    return pd.DataFrame(selected)[columns]


def fit_key(scenario: dict) -> tuple[int, str]:
    return int(scenario["training_count"]), str(scenario["kernel_family"])


def scenario_seed(config: dict, replicate: int, training_count: int, kernel_family: str) -> int:
    kernel_offset = {"squared_exponential": 0, "matern_3_2": 1, "matern_5_2": 2}[kernel_family]
    return int(config["joint_posterior_sampling"]["base_random_seed"]) + 10000 * replicate + 100 * training_count + kernel_offset


def verify_frozen_inputs(config: dict):
    protocol_dir = ROOT / config["protocol_output_directory"]
    manifest_path = protocol_dir / "frozen_robustness_protocol_manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if not manifest.get("protocol_frozen") or not manifest.get("masked_replay"):
        raise ValueError("Robustness protocol is not frozen for masked replay")
    if sha256(CONFIG_FILE) != manifest["configuration_sha256"]:
        raise ValueError("Robustness configuration changed after protocol freeze")
    expected_script_hash = manifest["implementation_sha256"][Path(__file__).name]
    if sha256(Path(__file__)) != expected_script_hash:
        raise ValueError("Proposal implementation changed after protocol freeze")
    bank_path = ROOT / config["complete_label_bank_csv"]
    if sha256(bank_path) != manifest["complete_label_bank_sha256"]:
        raise ValueError("Complete label bank changed after protocol freeze")
    designs_path = protocol_dir / "nested_matched_initial_designs.csv"
    if sha256(designs_path) != manifest["nested_matched_initial_designs_sha256"]:
        raise ValueError("Nested matched initial designs changed after protocol freeze")
    return manifest_path, manifest, designs_path, bank_path


def main() -> None:
    started = time.perf_counter()
    config = json.loads(CONFIG_FILE.read_text(encoding="utf-8"))
    outdir = ROOT / config["proposal_output_directory"]
    if outdir.exists():
        raise FileExistsError(f"Refusing to overwrite frozen robustness proposals: {outdir}")
    protocol_manifest_path, protocol, designs_path, bank_path = verify_frozen_inputs(config)

    geometry = pd.read_csv(ROOT / config["geometry_manifest_csv"]).sort_values("state_index").reset_index(drop=True)
    designs = pd.read_csv(designs_path)
    # The bank is opened only to expose labels belonging to each current
    # training prefix. Candidate frames below are built from geometry alone.
    bank = pd.read_csv(bank_path, usecols=["state_index", "final_label"])
    model_cfg = config["model"]
    features = model_cfg["feature_columns"]
    scenarios = config["scenarios"]
    methods = config["methods"]
    unique_fit_keys = sorted({fit_key(scenario) for scenario in scenarios})

    proposal_frames = []
    training_frames = []
    hyperparameter_rows = []
    diagnostic_rows = []
    for replicate in range(1, int(config["matched_designs"]["replicate_count"]) + 1):
        full_design = designs[designs.replicate == replicate].sort_values("design_position")
        if len(full_design) != int(config["matched_designs"]["maximum_training_count"]):
            raise ValueError(f"Replicate {replicate} maximum design size changed")
        fit_cache = {}
        for training_count, kernel_family in unique_fit_keys:
            design = full_design.iloc[:training_count].copy()
            training = design.merge(bank, on="state_index", how="left", validate="one_to_one")
            if len(training) != training_count or set(training.final_label) != {"B", "C"}:
                raise ValueError(f"Replicate {replicate}, n0={training_count} lacks binary B/C training labels")
            training["training_count"] = training_count
            training["kernel_family"] = kernel_family
            training["label_exposed_during_fit"] = True
            training_frames.append(training)
            X_train = training[features].to_numpy(float)
            y_train = np.where(training.final_label == model_cfg["positive_class"], 1.0, -1.0)
            state, optimization = optimize_hyperparameters(
                X_train,
                y_train,
                model_cfg["initial_starts"],
                kernel_family=kernel_family,
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
                raise RuntimeError(f"GP Laplace fit failed for replicate {replicate}, n0={training_count}, {kernel_family}")
            candidate = geometry[~geometry.state_index.isin(training.state_index)].copy().reset_index(drop=True)
            if "final_label" in candidate.columns:
                raise RuntimeError("Candidate outcome leaked into the acquisition frame")
            pairs = candidate_pairs(candidate, config)
            common = common_pair_geometry(state, candidate, pairs, X_train, config)
            seed = scenario_seed(config, replicate, training_count, kernel_family)
            opposite, negative_count, minimum_eigenvalue = joint_opposite_probability(
                state,
                common["X"],
                common["first"],
                common["second"],
                config["joint_posterior_sampling"],
                seed,
            )
            fit_cache[(training_count, kernel_family)] = (candidate, pairs, common, opposite, X_train)
            hyperparameter_rows.append(
                {
                    "replicate": replicate,
                    "training_count": training_count,
                    "kernel_family": kernel_family,
                    "B_training_count": int((training.final_label == "B").sum()),
                    "C_training_count": int((training.final_label == "C").sum()),
                    "lengthscale_1": state.lengthscales[0],
                    "lengthscale_2": state.lengthscales[1],
                    "signal_std": state.signal_std,
                    "log_marginal_likelihood": state.log_marginal_likelihood,
                    "optimizer_start_count": len(model_cfg["initial_starts"]),
                    "optimizer_record_count": len(optimization),
                }
            )
            diagnostic_rows.append(
                {
                    "replicate": replicate,
                    "training_count": training_count,
                    "kernel_family": kernel_family,
                    "candidate_count": len(candidate),
                    "admissible_pair_count": len(pairs),
                    "joint_sample_seed": seed,
                    "negative_joint_covariance_eigenvalue_count_before_clipping": negative_count,
                    "minimum_joint_covariance_eigenvalue_before_clipping": minimum_eigenvalue,
                }
            )

        for scenario in scenarios:
            training_count, kernel_family = fit_key(scenario)
            _, pairs, common, opposite, _ = fit_cache[(training_count, kernel_family)]
            distance_scale = float(config["pair_geometry"]["baseline_distance_scale"]) * float(
                scenario["distance_scale_multiplier"]
            )
            for method in methods:
                ranking = score_ranking(
                    method,
                    str(scenario["pair_score"]),
                    distance_scale,
                    pairs,
                    common,
                    opposite,
                )
                batch = greedy_batch(ranking, config)
                batch.insert(0, "distance_scale", distance_scale)
                batch.insert(0, "pair_score", scenario["pair_score"])
                batch.insert(0, "kernel_family", kernel_family)
                batch.insert(0, "training_count", training_count)
                batch.insert(0, "method", method)
                batch.insert(0, "scenario_level", scenario["level"])
                batch.insert(0, "scenario_factor", scenario["factor"])
                batch.insert(0, "scenario_id", scenario["scenario_id"])
                batch.insert(0, "replicate", replicate)
                proposal_frames.append(batch)
        print(f"replicate {replicate:02d}/40: fit {len(unique_fit_keys)} GPs and froze {len(scenarios) * len(methods) * 4} pairs", flush=True)

    proposals = pd.concat(proposal_frames, ignore_index=True)
    training_labels = pd.concat(training_frames, ignore_index=True)
    expected = (
        int(config["matched_designs"]["replicate_count"])
        * len(scenarios)
        * len(methods)
        * int(config["batch"]["pair_count"])
    )
    if len(proposals) != expected:
        raise ValueError(f"Frozen proposal count {len(proposals)} differs from expected {expected}")
    endpoint_count = proposals.groupby(["replicate", "scenario_id", "method"]).apply(
        lambda frame: len(set(frame.first_state_index) | set(frame.second_state_index))
    )
    if not np.all(endpoint_count == int(config["batch"]["new_endpoint_label_budget"])):
        raise ValueError("A proposal batch does not use exactly eight disjoint endpoints")

    outdir.mkdir(parents=True)
    proposal_path = outdir / "frozen_robustness_pair_proposals.csv"
    training_path = outdir / "exposed_training_prefix_labels.csv"
    hyperparameter_path = outdir / "robustness_GP_hyperparameters.csv"
    diagnostics_path = outdir / "robustness_acquisition_numerical_diagnostics.csv"
    proposals.to_csv(proposal_path, index=False)
    training_labels.to_csv(training_path, index=False)
    pd.DataFrame(hyperparameter_rows).to_csv(hyperparameter_path, index=False)
    pd.DataFrame(diagnostic_rows).to_csv(diagnostics_path, index=False)
    manifest = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "proposal_batches_frozen": True,
        "masked_replay": True,
        "replicate_count": int(config["matched_designs"]["replicate_count"]),
        "scenario_count": len(scenarios),
        "methods": methods,
        "pair_count_per_batch": int(config["batch"]["pair_count"]),
        "endpoint_budget_per_batch": int(config["batch"]["new_endpoint_label_budget"]),
        "proposal_count": len(proposals),
        "unique_GP_fit_count": len(hyperparameter_rows),
        "candidate_outcomes_used_for_proposals": False,
        "only_current_training_prefix_labels_exposed": True,
        "exact_edge_orbit_used": False,
        "heldout_CBA_region_used": False,
        "new_PDE_integrations": 0,
        "frozen_robustness_pair_proposals_sha256": sha256(proposal_path),
        "exposed_training_prefix_labels_sha256": sha256(training_path),
        "robustness_GP_hyperparameters_sha256": sha256(hyperparameter_path),
        "robustness_acquisition_numerical_diagnostics_sha256": sha256(diagnostics_path),
        "complete_label_bank_sha256": protocol["complete_label_bank_sha256"],
        "frozen_robustness_protocol_manifest_sha256": sha256(protocol_manifest_path),
        "proposal_script_sha256": sha256(Path(__file__)),
        "runtime_seconds": time.perf_counter() - started,
        "scientific_scope": config["evaluation_scope"]["claim"],
    }
    manifest_path = outdir / "frozen_robustness_proposal_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(json.dumps(manifest, indent=2))
    print(f"Frozen robustness proposals saved in: {outdir}")


if __name__ == "__main__":
    main()
