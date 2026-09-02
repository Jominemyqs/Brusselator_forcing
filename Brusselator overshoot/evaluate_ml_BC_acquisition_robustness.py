#!/usr/bin/env python3
"""Open frozen replay outcomes and evaluate acquisition robustness."""

from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import binomtest


ROOT = Path(__file__).resolve().parent
CONFIG_FILE = ROOT / "ml_bc_acquisition_robustness_config.json"
METHOD_ORDER = ["space_filling", "pointwise_uncertainty", "pair_aware"]
METHOD_COLORS = {
    "space_filling": "#009E73",
    "pointwise_uncertainty": "#E69F00",
    "pair_aware": "#CC79A7",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def paired_bootstrap_interval(difference: np.ndarray, draws: int, seed: int) -> tuple[float, float]:
    difference = np.asarray(difference, dtype=float)
    rng = np.random.default_rng(seed)
    sample = rng.choice(difference, size=(draws, len(difference)), replace=True).mean(axis=1)
    return tuple(np.quantile(sample, [0.025, 0.975]))


def mean_bootstrap_interval(values: np.ndarray, draws: int, seed: int) -> tuple[float, float]:
    values = np.asarray(values, dtype=float)
    rng = np.random.default_rng(seed)
    sample = rng.choice(values, size=(draws, len(values)), replace=True).mean(axis=1)
    return tuple(np.quantile(sample, [0.025, 0.975]))


def wilson_interval(successes: int, total: int, z: float = 1.959963984540054) -> tuple[float, float]:
    if total <= 0:
        return np.nan, np.nan
    proportion = successes / total
    denominator = 1.0 + z**2 / total
    center = (proportion + z**2 / (2.0 * total)) / denominator
    half_width = z * np.sqrt(proportion * (1.0 - proportion) / total + z**2 / (4.0 * total**2)) / denominator
    return max(0.0, center - half_width), min(1.0, center + half_width)


def compare_paired(
    left: np.ndarray,
    right: np.ndarray,
    comparison: str,
    draws: int,
    seed: int,
) -> dict:
    difference = np.asarray(left, dtype=float) - np.asarray(right, dtype=float)
    lower, upper = paired_bootstrap_interval(difference, draws, seed)
    wins = int((difference > 0).sum())
    ties = int((difference == 0).sum())
    losses = int((difference < 0).sum())
    non_ties = wins + losses
    sign_pvalue = float(binomtest(wins, non_ties, 0.5, alternative="two-sided").pvalue) if non_ties else 1.0
    return {
        "comparison": comparison,
        "replicate_count": len(difference),
        "mean_paired_difference_in_valid_brackets": float(difference.mean()),
        "paired_bootstrap_95_lower": lower,
        "paired_bootstrap_95_upper": upper,
        "left_wins": wins,
        "ties": ties,
        "left_losses": losses,
        "exact_two_sided_sign_test_pvalue_ignoring_ties": sign_pvalue,
    }


def verify_frozen_inputs(config: dict):
    protocol_dir = ROOT / config["protocol_output_directory"]
    protocol_manifest_path = protocol_dir / "frozen_robustness_protocol_manifest.json"
    protocol = json.loads(protocol_manifest_path.read_text(encoding="utf-8"))
    proposal_dir = ROOT / config["proposal_output_directory"]
    proposal_manifest_path = proposal_dir / "frozen_robustness_proposal_manifest.json"
    proposal_manifest = json.loads(proposal_manifest_path.read_text(encoding="utf-8"))
    proposal_path = proposal_dir / "frozen_robustness_pair_proposals.csv"
    if not proposal_manifest.get("proposal_batches_frozen") or proposal_manifest.get("candidate_outcomes_used_for_proposals"):
        raise ValueError("Robustness proposals were not frozen under masked replay")
    if sha256(proposal_path) != proposal_manifest["frozen_robustness_pair_proposals_sha256"]:
        raise ValueError("Robustness proposal table changed before outcome opening")
    if sha256(protocol_manifest_path) != proposal_manifest["frozen_robustness_protocol_manifest_sha256"]:
        raise ValueError("Frozen protocol manifest changed before outcome opening")
    if sha256(Path(__file__)) != protocol["implementation_sha256"][Path(__file__).name]:
        raise ValueError("Evaluation implementation changed after protocol freeze")
    bank_path = ROOT / config["complete_label_bank_csv"]
    if sha256(bank_path) != protocol["complete_label_bank_sha256"]:
        raise ValueError("Complete label bank changed before outcome opening")
    return protocol_manifest_path, protocol, proposal_manifest_path, proposal_manifest, proposal_path, bank_path


def make_figure(
    path_png: Path,
    path_pdf: Path,
    replicate_metrics: pd.DataFrame,
    method_summary: pd.DataFrame,
    budget_efficiency: pd.DataFrame,
) -> None:
    figure, axes = plt.subplots(2, 2, figsize=(11.5, 8.0), constrained_layout=True)

    baseline = replicate_metrics[replicate_metrics.scenario_id == "baseline"]
    for position, method in enumerate(METHOD_ORDER, start=1):
        values = baseline.loc[baseline.method == method, "valid_BC_brackets"].to_numpy(float)
        jitter = np.linspace(-0.11, 0.11, len(values))
        axes[0, 0].scatter(position + jitter, values, color=METHOD_COLORS[method], s=20, alpha=0.75)
        axes[0, 0].plot([position - 0.19, position + 0.19], [values.mean(), values.mean()], color="k", linewidth=2)
    axes[0, 0].set_xticks(range(1, 4), ["space\nfilling", "pointwise\nuncertainty", "pair\naware"])
    axes[0, 0].set_ylabel("valid B--C brackets among four")
    axes[0, 0].set_ylim(-0.2, 4.2)
    axes[0, 0].set_title("(a) Forty matched reference designs")
    axes[0, 0].grid(axis="y", alpha=0.25)

    scenario_order = [
        "initial_n15",
        "baseline",
        "initial_n35",
        "kernel_matern32",
        "kernel_matern52",
        "score_q",
        "score_qD",
        "distance_0p75",
        "distance_1p25",
    ]
    labels = ["n=15", "reference", "n=35", "Matérn 3/2", "Matérn 5/2", "q", "qD", "0.75d", "1.25d"]
    pair_summary = method_summary[
        (method_summary.method == "pair_aware") & method_summary.scenario_id.isin(scenario_order)
    ].set_index("scenario_id").loc[scenario_order]
    y = np.arange(len(scenario_order))
    mean = pair_summary.mean_valid_brackets_per_replicate.to_numpy(float)
    lower = pair_summary.bootstrap_95_lower.to_numpy(float)
    upper = pair_summary.bootstrap_95_upper.to_numpy(float)
    axes[0, 1].errorbar(mean, y, xerr=np.vstack([mean - lower, upper - mean]), fmt="o", color=METHOD_COLORS["pair_aware"], capsize=3)
    axes[0, 1].set_yticks(y, labels)
    axes[0, 1].invert_yaxis()
    axes[0, 1].set_xlim(0, 4.05)
    axes[0, 1].set_xlabel("mean valid brackets among four")
    axes[0, 1].set_title("(b) Pair-aware sensitivity")
    axes[0, 1].grid(axis="x", alpha=0.25)

    budget = budget_efficiency[budget_efficiency.scenario_id == "baseline"]
    for method in METHOD_ORDER:
        frame = budget[budget.method == method].sort_values("revealed_endpoint_labels")
        axes[1, 0].plot(
            frame.revealed_endpoint_labels,
            frame.probability_at_least_one_valid_bracket,
            "-o",
            color=METHOD_COLORS[method],
            label=method.replace("_", " "),
        )
    axes[1, 0].set_xticks([2, 4, 6, 8])
    axes[1, 0].set_ylim(-0.03, 1.03)
    axes[1, 0].set_xlabel("revealed endpoint labels")
    axes[1, 0].set_ylabel("fraction of designs with a bracket")
    axes[1, 0].set_title("(c) Reference budget efficiency")
    axes[1, 0].legend(fontsize=8)
    axes[1, 0].grid(alpha=0.25)

    ablation_ids = ["score_q", "score_qD", "baseline"]
    ablation_labels = ["q", "qD", "qDO"]
    ablation = method_summary[
        (method_summary.method == "pair_aware") & method_summary.scenario_id.isin(ablation_ids)
    ].set_index("scenario_id").loc[ablation_ids]
    values = ablation.mean_valid_brackets_per_replicate.to_numpy(float)
    lower = ablation.bootstrap_95_lower.to_numpy(float)
    upper = ablation.bootstrap_95_upper.to_numpy(float)
    axes[1, 1].bar(np.arange(3), values, color=["#56B4E9", "#0072B2", METHOD_COLORS["pair_aware"]])
    axes[1, 1].errorbar(np.arange(3), values, yerr=np.vstack([values - lower, upper - values]), fmt="none", color="k", capsize=3)
    axes[1, 1].set_xticks(np.arange(3), ablation_labels)
    axes[1, 1].set_ylim(0, 4.05)
    axes[1, 1].set_ylabel("mean valid brackets among four")
    axes[1, 1].set_title("(d) Pair-score ablation")
    axes[1, 1].grid(axis="y", alpha=0.25)

    figure.suptitle("Masked-replay robustness of task-aware B--C bracket acquisition")
    figure.savefig(path_png, dpi=300)
    figure.savefig(path_pdf)
    plt.close(figure)


def main() -> None:
    config = json.loads(CONFIG_FILE.read_text(encoding="utf-8"))
    outdir = ROOT / config["evaluation_output_directory"]
    if outdir.exists():
        raise FileExistsError(f"Refusing to overwrite robustness evaluation: {outdir}")
    (
        protocol_manifest_path,
        protocol,
        proposal_manifest_path,
        proposal_manifest,
        proposal_path,
        bank_path,
    ) = verify_frozen_inputs(config)

    proposals = pd.read_csv(proposal_path)
    bank = pd.read_csv(bank_path, usecols=["state_index", "final_label", "raw_five_way_outcome"])
    first = bank.rename(
        columns={
            "state_index": "first_state_index",
            "final_label": "first_outcome",
            "raw_five_way_outcome": "first_raw_outcome",
        }
    )
    second = bank.rename(
        columns={
            "state_index": "second_state_index",
            "final_label": "second_outcome",
            "raw_five_way_outcome": "second_raw_outcome",
        }
    )
    opened = proposals.merge(first, on="first_state_index", validate="many_to_one").merge(
        second, on="second_state_index", validate="many_to_one"
    )
    opened["verified_BC_bracket"] = [
        {left, right} == {"B", "C"} for left, right in zip(opened.first_outcome, opened.second_outcome)
    ]
    opened["contains_other_or_U"] = ~opened.first_outcome.isin(["B", "C"]) | ~opened.second_outcome.isin(["B", "C"])
    if opened.contains_other_or_U.any():
        raise ValueError("The controlled replay plane is no longer binary B/C")

    replicate_rows = []
    for (replicate, scenario_id, method), frame in opened.groupby(["replicate", "scenario_id", "method"], sort=True):
        frame = frame.sort_values("proposal_rank")
        valid = frame.verified_BC_bracket.to_numpy(bool)
        valid_ranks = frame.loc[frame.verified_BC_bracket, "proposal_rank"]
        first = frame.iloc[0]
        replicate_rows.append(
            {
                "replicate": int(replicate),
                "scenario_id": scenario_id,
                "scenario_factor": first.scenario_factor,
                "scenario_level": first.scenario_level,
                "training_count": int(first.training_count),
                "kernel_family": first.kernel_family,
                "pair_score": first.pair_score,
                "distance_scale": float(first.distance_scale),
                "method": method,
                "proposal_count": len(frame),
                "new_endpoint_label_budget": len(set(frame.first_state_index) | set(frame.second_state_index)),
                "valid_BC_brackets": int(valid.sum()),
                "bracket_fraction": float(valid.mean()),
                "valid_brackets_per_endpoint_label": float(valid.sum() / int(config["batch"]["new_endpoint_label_budget"])),
                "first_valid_bracket_rank": float(valid_ranks.min()) if len(valid_ranks) else np.nan,
                "endpoint_labels_to_first_valid_bracket": float(2 * valid_ranks.min()) if len(valid_ranks) else np.nan,
            }
        )
    replicate_metrics = pd.DataFrame(replicate_rows)

    statistics = config["statistics"]
    draws = int(statistics["paired_bootstrap_draws"])
    base_seed = int(statistics["paired_bootstrap_seed"])
    method_rows = []
    grouped = replicate_metrics.groupby(["scenario_id", "method"], sort=False)
    for group_number, ((scenario_id, method), frame) in enumerate(grouped, start=1):
        first = frame.iloc[0]
        values = frame.valid_BC_brackets.to_numpy(float)
        lower, upper = mean_bootstrap_interval(values, draws, base_seed + group_number)
        method_rows.append(
            {
                "scenario_id": scenario_id,
                "scenario_factor": first.scenario_factor,
                "scenario_level": first.scenario_level,
                "training_count": int(first.training_count),
                "kernel_family": first.kernel_family,
                "pair_score": first.pair_score,
                "distance_scale": float(first.distance_scale),
                "method": method,
                "replicate_count": len(frame),
                "proposal_count": int(frame.proposal_count.sum()),
                "algorithmic_endpoint_label_budget": int(frame.new_endpoint_label_budget.sum()),
                "total_valid_BC_brackets": int(frame.valid_BC_brackets.sum()),
                "pooled_bracket_yield": float(frame.valid_BC_brackets.sum() / frame.proposal_count.sum()),
                "mean_valid_brackets_per_replicate": float(values.mean()),
                "bootstrap_95_lower": lower,
                "bootstrap_95_upper": upper,
                "replicates_with_at_least_one_bracket": int((values > 0).sum()),
                "coverage_fraction": float((values > 0).mean()),
                "valid_brackets_per_endpoint_label": float(frame.valid_BC_brackets.sum() / frame.new_endpoint_label_budget.sum()),
            }
        )
    method_summary = pd.DataFrame(method_rows)

    comparison_rows = []
    for scenario_number, scenario in enumerate(config["scenarios"]):
        scenario_id = scenario["scenario_id"]
        frame = replicate_metrics[replicate_metrics.scenario_id == scenario_id]
        wide = frame.pivot(index="replicate", columns="method", values="valid_BC_brackets")
        for comparator_number, comparator in enumerate(("pointwise_uncertainty", "space_filling")):
            result = compare_paired(
                wide.pair_aware.to_numpy(float),
                wide[comparator].to_numpy(float),
                f"pair_aware_minus_{comparator}",
                draws,
                base_seed + 1000 + 10 * scenario_number + comparator_number,
            )
            result.update({"scenario_id": scenario_id, "scenario_factor": scenario["factor"], "scenario_level": scenario["level"]})
            comparison_rows.append(result)
    paired_comparisons = pd.DataFrame(comparison_rows)

    ablation_map = {"q": "score_q", "qD": "score_qD", "qDO": "baseline"}
    ablation = replicate_metrics[
        (replicate_metrics.method == "pair_aware") & replicate_metrics.scenario_id.isin(ablation_map.values())
    ].copy()
    inverse = {value: key for key, value in ablation_map.items()}
    ablation["ablation_score"] = ablation.scenario_id.map(inverse)
    ablation_wide = ablation.pivot(index="replicate", columns="ablation_score", values="valid_BC_brackets")
    ablation_rows = []
    for index, (left, right) in enumerate((("qDO", "qD"), ("qD", "q"), ("qDO", "q"))):
        result = compare_paired(
            ablation_wide[left].to_numpy(float),
            ablation_wide[right].to_numpy(float),
            f"{left}_minus_{right}",
            draws,
            base_seed + 2000 + index,
        )
        ablation_rows.append(result)
    ablation_comparisons = pd.DataFrame(ablation_rows)

    budget_rows = []
    for (scenario_id, method), frame in opened.groupby(["scenario_id", "method"], sort=False):
        for endpoint_labels in statistics["budget_curve_endpoint_labels"]:
            maximum_rank = int(endpoint_labels) // 2
            success_by_replicate = frame[frame.proposal_rank <= maximum_rank].groupby("replicate").verified_BC_bracket.any()
            successes = int(success_by_replicate.sum())
            total = len(success_by_replicate)
            lower, upper = wilson_interval(successes, total)
            budget_rows.append(
                {
                    "scenario_id": scenario_id,
                    "method": method,
                    "revealed_endpoint_labels": int(endpoint_labels),
                    "successful_replicates": successes,
                    "replicate_count": total,
                    "probability_at_least_one_valid_bracket": successes / total,
                    "wilson_95_lower": lower,
                    "wilson_95_upper": upper,
                }
            )
    budget_efficiency = pd.DataFrame(budget_rows)

    outdir.mkdir(parents=True)
    paths = {
        "opened_proposals": outdir / "robustness_proposals_with_opened_outcomes.csv",
        "replicate_metrics": outdir / "robustness_replicate_metrics.csv",
        "method_summary": outdir / "robustness_method_summary.csv",
        "paired_comparisons": outdir / "robustness_paired_comparisons.csv",
        "ablation_comparisons": outdir / "robustness_ablation_comparisons.csv",
        "budget_efficiency": outdir / "robustness_budget_efficiency.csv",
    }
    opened.to_csv(paths["opened_proposals"], index=False)
    replicate_metrics.to_csv(paths["replicate_metrics"], index=False)
    method_summary.to_csv(paths["method_summary"], index=False)
    paired_comparisons.to_csv(paths["paired_comparisons"], index=False)
    ablation_comparisons.to_csv(paths["ablation_comparisons"], index=False)
    budget_efficiency.to_csv(paths["budget_efficiency"], index=False)
    figure_png = outdir / "ml_BC_acquisition_robustness_overview.png"
    figure_pdf = outdir / "ml_BC_acquisition_robustness_overview.pdf"
    make_figure(figure_png, figure_pdf, replicate_metrics, method_summary, budget_efficiency)

    baseline_results = method_summary[method_summary.scenario_id == "baseline"].to_dict(orient="records")
    summary = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "status": "complete masked-replay robustness evaluation",
        "new_PDE_integrations": 0,
        "replicate_count": int(config["matched_designs"]["replicate_count"]),
        "scenario_count": len(config["scenarios"]),
        "primary_independent_unit": statistics["primary_unit"],
        "candidate_outcomes_opened_only_after_proposal_hash_freeze": True,
        "baseline_results": baseline_results,
        "baseline_paired_results": paired_comparisons[paired_comparisons.scenario_id == "baseline"].to_dict(orient="records"),
        "pair_score_ablation_results": method_summary[
            (method_summary.method == "pair_aware") & method_summary.scenario_id.isin(ablation_map.values())
        ].to_dict(orient="records"),
        "pair_score_ablation_paired_comparisons": ablation_comparisons.to_dict(orient="records"),
        "integrity": {
            "frozen_robustness_protocol_manifest_sha256": sha256(protocol_manifest_path),
            "frozen_robustness_proposal_manifest_sha256": sha256(proposal_manifest_path),
            "complete_label_bank_sha256": protocol["complete_label_bank_sha256"],
            **{f"{name}_sha256": sha256(path) for name, path in paths.items()},
            "figure_png_sha256": sha256(figure_png),
            "figure_pdf_sha256": sha256(figure_pdf),
        },
        "scientific_scope": config["evaluation_scope"]["claim"],
        "excluded_claims": config["evaluation_scope"]["excluded_claims"],
    }
    summary_path = outdir / "ml_BC_acquisition_robustness_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(method_summary.to_string(index=False))
    print(paired_comparisons.to_string(index=False))
    print(ablation_comparisons.to_string(index=False))
    print(f"Robustness evaluation saved in: {outdir}")


if __name__ == "__main__":
    main()
