#!/usr/bin/env python3
"""Post hoc orientation ablation for the 40-design pointwise comparator.

This masked replay uses the already frozen development geometry, nested
initial designs, GP specification, and complete PDE label bank.  It changes
only the pointwise score from sqrt(H_i H_j) D_ij O_ij to
sqrt(H_i H_j) D_ij.  The calculation is an algorithmic comparator audit; it
is not prospective evidence and performs no new PDE integrations.
"""

from __future__ import annotations

import hashlib
import json
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import binomtest

from probit_gp_kernels import optimize_hyperparameters
from run_ml_BC_acquisition_robustness_proposals import (
    candidate_pairs,
    common_pair_geometry,
    greedy_batch,
)


ROOT = Path(__file__).resolve().parent
CONFIG_FILE = ROOT / "ml_bc_acquisition_qd_primary_config.json"
OUTPUT_DIRECTORY = ROOT / "experiment_outputs/ml_BC_pointwise_orientation_ablation_v1"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def pointwise_orientation_free_ranking(
    pairs: pd.DataFrame,
    common: dict,
    distance_scale: float,
) -> pd.DataFrame:
    first = common["first"]
    second = common["second"]
    distance_penalty = np.exp(-0.5 * (common["separation"] / distance_scale) ** 2)
    information = np.sqrt(common["entropy"][first] * common["entropy"][second])
    ranked = pairs.copy()
    ranked["first_probability_C"] = common["probability"][first]
    ranked["second_probability_C"] = common["probability"][second]
    ranked["first_entropy_bits"] = common["entropy"][first]
    ranked["second_entropy_bits"] = common["entropy"][second]
    ranked["first_novelty"] = common["novelty"][first]
    ranked["second_novelty"] = common["novelty"][second]
    ranked["distance_penalty"] = distance_penalty
    ranked["raw_normal_orientation"] = common["raw_orientation"]
    ranked["orientation_factor"] = common["orientation"]
    ranked["midpoint_gradient_norm"] = common["gradient_norm"]
    ranked["information_factor"] = information
    ranked["acquisition_score"] = information * distance_penalty
    ranked["score_definition"] = "sqrt(endpoint_entropy_product)*D"
    ranked = ranked.sort_values(
        ["acquisition_score", "information_factor", "normalized_separation", "pair_id"],
        ascending=[False, False, True, True],
    ).reset_index(drop=True)
    ranked.insert(1, "score_rank", np.arange(1, len(ranked) + 1))
    return ranked


def bootstrap_mean_interval(values: np.ndarray, seed: int, draws: int = 20_000) -> tuple[float, float]:
    rng = np.random.default_rng(seed)
    samples = rng.choice(values, size=(draws, len(values)), replace=True).mean(axis=1)
    return tuple(np.quantile(samples, [0.025, 0.975]))


def paired_comparison(left: np.ndarray, right: np.ndarray, seed: int) -> dict:
    difference = np.asarray(left, float) - np.asarray(right, float)
    lower, upper = bootstrap_mean_interval(difference, seed)
    wins = int((difference > 0).sum())
    ties = int((difference == 0).sum())
    losses = int((difference < 0).sum())
    non_ties = wins + losses
    pvalue = float(binomtest(wins, non_ties, 0.5).pvalue) if non_ties else 1.0
    return {
        "mean_paired_difference": float(difference.mean()),
        "paired_bootstrap_95_lower": lower,
        "paired_bootstrap_95_upper": upper,
        "left_wins": wins,
        "ties": ties,
        "left_losses": losses,
        "exact_two_sided_sign_test_pvalue_ignoring_ties": pvalue,
    }


def main() -> None:
    started = time.perf_counter()
    if OUTPUT_DIRECTORY.exists():
        raise FileExistsError(f"Refusing to overwrite {OUTPUT_DIRECTORY}")
    config = json.loads(CONFIG_FILE.read_text(encoding="utf-8"))
    protocol_dir = ROOT / config["protocol_output_directory"]
    protocol_manifest_path = protocol_dir / "frozen_robustness_protocol_manifest.json"
    protocol = json.loads(protocol_manifest_path.read_text(encoding="utf-8"))
    designs_path = protocol_dir / "nested_matched_initial_designs.csv"
    geometry_path = ROOT / config["geometry_manifest_csv"]
    bank_path = ROOT / config["complete_label_bank_csv"]
    if sha256(bank_path) != protocol["complete_label_bank_sha256"]:
        raise ValueError("Complete label bank differs from the frozen qD-primary protocol")
    if sha256(designs_path) != protocol["nested_matched_initial_designs_sha256"]:
        raise ValueError("Matched initial designs differ from the frozen qD-primary protocol")

    geometry = pd.read_csv(geometry_path).sort_values("state_index").reset_index(drop=True)
    designs = pd.read_csv(designs_path)
    bank = pd.read_csv(bank_path, usecols=["state_index", "final_label"])
    model_cfg = config["model"]
    features = model_cfg["feature_columns"]
    training_count = 25
    distance_scale = float(config["pair_geometry"]["baseline_distance_scale"])

    proposal_frames = []
    fit_rows = []
    for replicate in range(1, int(config["matched_designs"]["replicate_count"]) + 1):
        design = (
            designs[designs.replicate == replicate]
            .sort_values("design_position")
            .iloc[:training_count]
            .merge(bank, on="state_index", validate="one_to_one")
        )
        if len(design) != training_count or set(design.final_label) != {"B", "C"}:
            raise ValueError(f"Replicate {replicate} lacks the declared binary training design")
        X_train = design[features].to_numpy(float)
        y_train = np.where(design.final_label == model_cfg["positive_class"], 1.0, -1.0)
        state, optimization = optimize_hyperparameters(
            X_train,
            y_train,
            model_cfg["initial_starts"],
            kernel_family="squared_exponential",
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
            raise RuntimeError(f"GP fit failed for replicate {replicate}")
        candidate = geometry[~geometry.state_index.isin(design.state_index)].copy().reset_index(drop=True)
        pairs = candidate_pairs(candidate, config)
        common = common_pair_geometry(state, candidate, pairs, X_train, config)
        ranking = pointwise_orientation_free_ranking(pairs, common, distance_scale)
        batch = greedy_batch(ranking, config)
        batch.insert(0, "method", "pointwise_uncertainty_no_orientation")
        batch.insert(0, "replicate", replicate)
        proposal_frames.append(batch)
        fit_rows.append(
            {
                "replicate": replicate,
                "lengthscale_1": state.lengthscales[0],
                "lengthscale_2": state.lengthscales[1],
                "signal_std": state.signal_std,
                "optimizer_record_count": len(optimization),
            }
        )
        print(f"replicate {replicate:02d}/40 complete", flush=True)

    proposals = pd.concat(proposal_frames, ignore_index=True)
    first = bank.rename(columns={"state_index": "first_state_index", "final_label": "first_outcome"})
    second = bank.rename(columns={"state_index": "second_state_index", "final_label": "second_outcome"})
    opened = proposals.merge(first, on="first_state_index", validate="many_to_one").merge(
        second, on="second_state_index", validate="many_to_one"
    )
    opened["opposite_label_BC_pair"] = [
        {left, right} == {"B", "C"} for left, right in zip(opened.first_outcome, opened.second_outcome)
    ]
    replicate = (
        opened.groupby("replicate", as_index=False)
        .agg(
            opposite_label_pair_count=("opposite_label_BC_pair", "sum"),
            mean_pair_separation=("normalized_separation", "mean"),
            median_pair_separation=("normalized_separation", "median"),
        )
    )
    values = replicate.opposite_label_pair_count.to_numpy(float)
    lower, upper = bootstrap_mean_interval(values, 20261001)

    existing_root = ROOT / config["evaluation_output_directory"]
    existing_replicates = pd.read_csv(existing_root / "qD_primary_replicate_metrics.csv")
    baseline = existing_replicates[existing_replicates.scenario_id == "baseline"].pivot(
        index="replicate", columns="method", values="valid_BC_brackets"
    )
    merged = replicate.set_index("replicate").join(baseline, how="inner")
    comparisons = {
        "pair_aware_qD_minus_pointwise_no_orientation": paired_comparison(
            merged.pair_aware.to_numpy(float), merged.opposite_label_pair_count.to_numpy(float), 20261002
        ),
        "pointwise_with_orientation_minus_without_orientation": paired_comparison(
            merged.pointwise_uncertainty.to_numpy(float),
            merged.opposite_label_pair_count.to_numpy(float),
            20261003,
        ),
    }
    summary = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "status": "complete post hoc masked-replay pointwise orientation ablation",
        "new_PDE_integrations": 0,
        "replicate_count": int(len(replicate)),
        "proposal_count": int(len(opened)),
        "total_opposite_label_BC_pairs": int(values.sum()),
        "mean_pairs_per_replicate": float(values.mean()),
        "bootstrap_95_interval_for_mean": [lower, upper],
        "replicates_with_at_least_one_pair": int((values > 0).sum()),
        "coverage_fraction": float((values > 0).mean()),
        "pairs_per_revealed_endpoint_label": float(values.sum() / (8 * len(replicate))),
        "mean_proposed_separation": float(opened.normalized_separation.mean()),
        "median_proposed_separation": float(opened.normalized_separation.median()),
        "comparisons": comparisons,
        "interpretation": (
            "post hoc orientation-factor ablation on the same prelocalized development family; "
            "not prospective transfer evidence or a wall-clock savings measurement"
        ),
        "source_hashes": {
            "configuration": sha256(CONFIG_FILE),
            "protocol_manifest": sha256(protocol_manifest_path),
            "nested_designs": sha256(designs_path),
            "geometry_manifest": sha256(geometry_path),
            "complete_label_bank": sha256(bank_path),
            "implementation": sha256(Path(__file__)),
        },
        "runtime_seconds": time.perf_counter() - started,
    }

    OUTPUT_DIRECTORY.mkdir(parents=True)
    proposals.to_csv(OUTPUT_DIRECTORY / "pointwise_no_orientation_proposals.csv", index=False)
    opened.to_csv(OUTPUT_DIRECTORY / "pointwise_no_orientation_opened_outcomes.csv", index=False)
    replicate.to_csv(OUTPUT_DIRECTORY / "pointwise_no_orientation_replicate_metrics.csv", index=False)
    pd.DataFrame(fit_rows).to_csv(OUTPUT_DIRECTORY / "pointwise_no_orientation_gp_fits.csv", index=False)
    (OUTPUT_DIRECTORY / "pointwise_orientation_ablation_summary.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8"
    )
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
