#!/usr/bin/env python3
"""Fit the transfer GP and freeze four pair proposals per method."""

from __future__ import annotations

import hashlib
import json
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist

from probit_gp_kernels import optimize_hyperparameters
from run_ml_BC_acquisition_robustness_proposals import (
    candidate_pairs,
    common_pair_geometry,
    greedy_batch,
    joint_opposite_probability,
)


ROOT = Path(__file__).resolve().parent
CONFIG = ROOT / "ml_bc_transfer_prospective_config.json"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def entropy_bits(probability: np.ndarray) -> np.ndarray:
    probability = np.clip(probability, 1e-12, 1.0 - 1e-12)
    return -(probability * np.log2(probability) + (1.0 - probability) * np.log2(1.0 - probability))


def rank_method(method: str, pairs: pd.DataFrame, common: dict, opposite: np.ndarray, config: dict) -> pd.DataFrame:
    first = common["first"]
    second = common["second"]
    separation = common["separation"]
    distance_scale = float(config["pair_geometry"]["distance_scale"])
    distance_penalty = np.exp(-0.5 * (separation / distance_scale) ** 2)
    ranked = pairs.copy()
    ranked["first_probability_C"] = common["probability"][first]
    ranked["second_probability_C"] = common["probability"][second]
    ranked["first_entropy_bits"] = common["entropy"][first]
    ranked["second_entropy_bits"] = common["entropy"][second]
    ranked["first_novelty"] = common["novelty"][first]
    ranked["second_novelty"] = common["novelty"][second]
    ranked["joint_opposite_sign_probability"] = opposite
    ranked["distance_penalty"] = distance_penalty
    ranked["orientation_factor"] = common["orientation"]
    if method == "space_filling":
        information = np.sqrt(common["novelty"][first] * common["novelty"][second])
        definition = "sqrt(endpoint_novelty_product)*D"
    elif method == "pointwise_uncertainty":
        information = np.sqrt(common["entropy"][first] * common["entropy"][second])
        definition = "sqrt(endpoint_entropy_product)*D"
    elif method == "pair_aware":
        information = opposite
        definition = "qD"
    else:
        raise ValueError(f"Unknown method: {method}")
    ranked["information_factor"] = information
    ranked["acquisition_score"] = information * distance_penalty
    ranked["score_definition"] = definition
    return ranked.sort_values(
        ["acquisition_score", "information_factor", "normalized_separation", "pair_id"],
        ascending=[False, False, True, True],
    ).reset_index(drop=True)


def main() -> None:
    started = time.perf_counter()
    config = json.loads(CONFIG.read_text(encoding="utf-8"))
    outdir = ROOT / config["proposal_directory"]
    if outdir.exists():
        raise FileExistsError(f"Refusing to overwrite frozen transfer proposals: {outdir}")
    protocol_path = ROOT / config["protocol_directory"] / "frozen_transfer_protocol_manifest.json"
    protocol = json.loads(protocol_path.read_text(encoding="utf-8"))
    geometry_path = ROOT / config["geometry_manifest_csv"]
    initial_path = ROOT / config["initial_labels_csv"]
    if not protocol.get("protocol_frozen") or protocol.get("candidate_outcomes_used"):
        raise ValueError("Transfer protocol is not frozen")
    if sha256(geometry_path) != protocol["geometry_manifest_sha256"]:
        raise ValueError("Transfer geometry changed after protocol freeze")

    geometry = pd.read_csv(geometry_path).sort_values("state_index").reset_index(drop=True)
    initial = pd.read_csv(initial_path)
    if len(initial) != config["geometry_gate"]["expected_initial_count"]:
        raise ValueError("Initial transfer label count changed")
    if set(initial.final_label) != {"B", "C"} or not initial.final_label.isin(["B", "C"]).all():
        raise RuntimeError("Prospective transfer family failed the binary B/C initial-design gate")
    design_indices = set(initial.state_index.astype(int))
    expected_design = set(geometry.loc[geometry.design_role == "initial_seed", "state_index"].astype(int))
    if design_indices != expected_design:
        raise ValueError("Initial labels do not match the frozen design")

    model_cfg = config["model"]
    features = model_cfg["feature_columns"]
    training = geometry.merge(initial[["state_index", "final_label"]], on="state_index", how="inner", validate="one_to_one")
    X_train = training[features].to_numpy(float)
    y_train = np.where(training.final_label == model_cfg["positive_class"], 1.0, -1.0)
    state, optimization = optimize_hyperparameters(
        X_train,
        y_train,
        model_cfg["initial_starts"],
        kernel_family=model_cfg["kernel_family"],
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
        raise RuntimeError("Transfer GP did not converge")

    candidate = geometry[~geometry.state_index.isin(design_indices)].copy().reset_index(drop=True)
    candidate = candidate.drop(columns=[column for column in ["outcome_label"] if column in candidate], errors="ignore")
    pairs = candidate_pairs(candidate, config)
    common = common_pair_geometry(state, candidate, pairs, X_train, config)
    opposite, negative_count, minimum_eigenvalue = joint_opposite_probability(
        state,
        common["X"],
        common["first"],
        common["second"],
        config["joint_posterior_sampling"],
        int(config["joint_posterior_sampling"]["random_seed"]),
    )

    batches = []
    for method in config["methods"]:
        ranking = rank_method(method, pairs, common, opposite, config)
        batch = greedy_batch(ranking, config)
        batch.insert(0, "method", method)
        batches.append(batch)
    proposals = pd.concat(batches, ignore_index=True)
    if not np.all(proposals.groupby("method").size() == config["batch"]["pair_count"]):
        raise ValueError("A method did not produce exactly four pairs")
    endpoint_rows = []
    for row in proposals.itertuples(index=False):
        for side in ("first", "second"):
            endpoint_rows.append(
                {
                    "state_index": int(getattr(row, f"{side}_state_index")),
                    "sample_id": getattr(row, f"{side}_sample_id"),
                    "selected_by_method": row.method,
                    "proposal_rank": int(row.proposal_rank),
                    "pair_id": row.pair_id,
                    "pair_side": side,
                }
            )
    endpoint_usage = pd.DataFrame(endpoint_rows)
    endpoint_manifest = (
        endpoint_usage.groupby(["state_index", "sample_id"], as_index=False)
        .agg(
            selecting_methods=("selected_by_method", lambda values: ";".join(sorted(set(values)))),
            selection_count=("selected_by_method", "size"),
        )
        .sort_values("state_index")
    )

    outdir.mkdir(parents=True)
    proposals_path = outdir / "frozen_transfer_pair_proposals.csv"
    endpoints_path = outdir / "frozen_transfer_endpoint_manifest.csv"
    proposals.to_csv(proposals_path, index=False)
    endpoint_manifest.to_csv(endpoints_path, index=False)
    endpoint_usage.to_csv(outdir / "frozen_transfer_endpoint_usage.csv", index=False)
    pd.DataFrame(
        [
            {
                "lengthscale_1": state.lengthscales[0],
                "lengthscale_2": state.lengthscales[1],
                "signal_std": state.signal_std,
                "log_marginal_likelihood": state.log_marginal_likelihood,
                "optimizer_start_count": len(optimization),
                "negative_joint_covariance_eigenvalues_before_clipping": negative_count,
                "minimum_joint_covariance_eigenvalue_before_clipping": minimum_eigenvalue,
            }
        ]
    ).to_csv(outdir / "transfer_gp_diagnostics.csv", index=False)
    manifest = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "proposal_batches_frozen": True,
        "protocol_manifest_sha256": sha256(protocol_path),
        "geometry_manifest_sha256": sha256(geometry_path),
        "initial_labels_sha256": sha256(initial_path),
        "configuration_sha256": sha256(CONFIG),
        "proposal_table_sha256": sha256(proposals_path),
        "endpoint_manifest_sha256": sha256(endpoints_path),
        "implementation_sha256": {Path(__file__).name: sha256(Path(__file__))},
        "initial_labels_used_for_fit": True,
        "candidate_outcomes_used_for_proposals": False,
        "endpoint_outcomes_used_for_proposals": False,
        "exact_edge_orbit_used": False,
        "development_label_bank_used": False,
        "pair_count_per_method": config["batch"]["pair_count"],
        "endpoint_budget_per_method": config["batch"]["endpoint_budget_per_method"],
        "unique_selected_endpoint_count": len(endpoint_manifest),
        "methods": config["methods"],
        "pair_aware_score": "qD",
        "elapsed_algorithm_seconds": time.perf_counter() - started,
        "phase": "prospective selected-endpoint PDE labeling on an independent state family",
    }
    manifest_path = outdir / "frozen_transfer_proposal_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(proposals[["method", "proposal_rank", "pair_id", "normalized_separation", "acquisition_score"]].to_string(index=False))
    print(json.dumps(manifest, indent=2))
    print(f"Frozen transfer proposals saved in: {outdir}")


if __name__ == "__main__":
    main()
