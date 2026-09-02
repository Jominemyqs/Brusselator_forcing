#!/usr/bin/env python3
"""Open frozen endpoint outcomes and evaluate the matched acquisition study."""

from __future__ import annotations

import hashlib
import itertools
import json
from datetime import datetime, timezone
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parent
CONFIG_FILE = ROOT / "ml_bc_matched_benchmark_protocol_config.json"
BANK_DIR = ROOT / "experiment_outputs" / "ml_BC_matched_label_bank_v1"
PROPOSAL_DIR = ROOT / "experiment_outputs" / "ml_BC_matched_acquisition_proposals_v1"
OUTDIR = ROOT / "experiment_outputs" / "ml_BC_matched_acquisition_evaluation_v1"
METHOD_ORDER = ["space_filling", "pointwise_uncertainty", "pair_aware"]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def paired_sign_flip_pvalue(difference: np.ndarray) -> float:
    difference = np.asarray(difference, dtype=float)
    observed = abs(float(difference.mean()))
    values = []
    for signs in itertools.product((-1.0, 1.0), repeat=len(difference)):
        values.append(abs(float(np.mean(difference * np.asarray(signs)))))
    return float(np.mean(np.asarray(values) >= observed - 1e-15))


def paired_bootstrap_interval(difference: np.ndarray, seed: int = 20260837) -> tuple[float, float]:
    rng = np.random.default_rng(seed)
    difference = np.asarray(difference, dtype=float)
    draws = rng.choice(difference, size=(20000, len(difference)), replace=True).mean(axis=1)
    return tuple(np.quantile(draws, [0.025, 0.975]))


def make_figure(path: Path, replicate_metrics: pd.DataFrame, method_summary: pd.DataFrame) -> None:
    colors = {"space_filling": "#009E73", "pointwise_uncertainty": "#E69F00", "pair_aware": "#CC79A7"}
    figure, axes = plt.subplots(1, 3, figsize=(12, 4), constrained_layout=True)
    for position, method in enumerate(METHOD_ORDER, start=1):
        values = replicate_metrics.loc[replicate_metrics.method == method, "valid_BC_brackets"].to_numpy()
        jitter = np.linspace(-0.08, 0.08, len(values))
        axes[0].scatter(position + jitter, values, color=colors[method], s=28)
        axes[0].plot([position - 0.18, position + 0.18], [values.mean(), values.mean()], color="k", linewidth=2)
    axes[0].set_xticks(range(1, 4), ["space\nfilling", "pointwise\nuncertainty", "pair\naware"])
    axes[0].set_ylabel("valid B-C brackets among four")
    axes[0].set_ylim(-0.2, 4.2)
    axes[0].grid(axis="y", alpha=0.25)

    axes[1].bar(
        range(3),
        method_summary.total_valid_BC_brackets,
        color=[colors[item] for item in METHOD_ORDER],
    )
    axes[1].set_xticks(range(3), ["space\nfilling", "pointwise\nuncertainty", "pair\naware"])
    axes[1].set_ylabel("valid brackets among 48 proposals")
    axes[1].grid(axis="y", alpha=0.25)

    for method in METHOD_ORDER:
        frame = replicate_metrics[replicate_metrics.method == method]
        axes[2].plot(
            frame.replicate,
            frame.valid_BC_brackets,
            "-o",
            label=method.replace("_", " "),
            color=colors[method],
            markersize=4,
        )
    axes[2].set_xlabel("matched initial-design replicate")
    axes[2].set_ylabel("valid B-C brackets")
    axes[2].set_xticks(range(1, 13))
    axes[2].set_ylim(-0.2, 4.2)
    axes[2].legend(fontsize=8)
    axes[2].grid(alpha=0.25)
    figure.suptitle("Matched outcome-blind acquisition benchmark")
    figure.savefig(path, dpi=300)
    plt.close(figure)


def main() -> None:
    if OUTDIR.exists():
        raise FileExistsError(f"Refusing to overwrite matched evaluation: {OUTDIR}")
    config = json.loads(CONFIG_FILE.read_text(encoding="utf-8"))
    proposal_manifest_path = PROPOSAL_DIR / "frozen_proposal_manifest.json"
    proposal_manifest = json.loads(proposal_manifest_path.read_text(encoding="utf-8"))
    proposal_path = PROPOSAL_DIR / "frozen_matched_pair_proposals.csv"
    if not proposal_manifest.get("proposal_batches_frozen") or proposal_manifest.get("candidate_outcomes_used"):
        raise ValueError("Proposal batches were not frozen outcome-blind")
    if sha256(proposal_path) != proposal_manifest["frozen_matched_pair_proposals_sha256"]:
        raise ValueError("Frozen proposal table changed before evaluation")
    bank_manifest_path = BANK_DIR / "complete_label_bank_manifest.json"
    bank_manifest = json.loads(bank_manifest_path.read_text(encoding="utf-8"))
    bank_path = BANK_DIR / "complete_label_bank.csv"
    if sha256(bank_path) != bank_manifest["complete_label_bank_sha256"]:
        raise ValueError("Complete label bank changed before outcome reveal")

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
    opened["same_outcome_pair"] = opened.first_outcome == opened.second_outcome
    opened["contains_other_or_U"] = ~opened.first_outcome.isin(["B", "C"]) | ~opened.second_outcome.isin(
        ["B", "C"]
    )

    replicate_rows = []
    for (replicate, method), frame in opened.groupby(["replicate", "method"], sort=True):
        valid = frame.verified_BC_bracket.to_numpy(bool)
        ranks = frame.loc[frame.verified_BC_bracket, "proposal_rank"]
        replicate_rows.append(
            {
                "replicate": int(replicate),
                "method": method,
                "proposal_count": len(frame),
                "new_endpoint_label_budget": len(set(frame.first_state_index) | set(frame.second_state_index)),
                "valid_BC_brackets": int(valid.sum()),
                "bracket_fraction": float(valid.mean()),
                "valid_brackets_per_endpoint_label": float(valid.sum() / 8.0),
                "first_valid_bracket_rank": float(ranks.min()) if len(ranks) else np.nan,
                "endpoint_labels_to_first_valid_bracket": float(2 * ranks.min()) if len(ranks) else np.nan,
                "other_or_U_pair_count": int(frame.contains_other_or_U.sum()),
            }
        )
    replicate_metrics = pd.DataFrame(replicate_rows)

    method_rows = []
    for method in METHOD_ORDER:
        frame = replicate_metrics[replicate_metrics.method == method]
        method_rows.append(
            {
                "method": method,
                "replicate_count": len(frame),
                "proposal_count": int(frame.proposal_count.sum()),
                "algorithmic_endpoint_label_budget": int(frame.new_endpoint_label_budget.sum()),
                "total_valid_BC_brackets": int(frame.valid_BC_brackets.sum()),
                "pooled_bracket_fraction": float(frame.valid_BC_brackets.sum() / frame.proposal_count.sum()),
                "mean_valid_brackets_per_replicate": float(frame.valid_BC_brackets.mean()),
                "standard_deviation_across_replicates": float(frame.valid_BC_brackets.std(ddof=1)),
                "median_valid_brackets_per_replicate": float(frame.valid_BC_brackets.median()),
                "replicates_with_at_least_one_bracket": int((frame.valid_BC_brackets > 0).sum()),
                "mean_valid_brackets_per_endpoint_label": float(frame.valid_brackets_per_endpoint_label.mean()),
                "median_endpoint_labels_to_first_bracket_when_found": float(
                    frame.endpoint_labels_to_first_valid_bracket.dropna().median()
                )
                if frame.endpoint_labels_to_first_valid_bracket.notna().any()
                else np.nan,
            }
        )
    method_summary = pd.DataFrame(method_rows)

    wide = replicate_metrics.pivot(index="replicate", columns="method", values="valid_BC_brackets")
    comparison_rows = []
    for baseline in ("space_filling", "pointwise_uncertainty"):
        difference = wide.pair_aware.to_numpy(float) - wide[baseline].to_numpy(float)
        lower, upper = paired_bootstrap_interval(difference, seed=20260837 + len(comparison_rows))
        comparison_rows.append(
            {
                "comparison": f"pair_aware_minus_{baseline}",
                "replicate_count": len(difference),
                "mean_paired_difference_in_valid_brackets": float(difference.mean()),
                "paired_bootstrap_95_lower": lower,
                "paired_bootstrap_95_upper": upper,
                "exact_sign_flip_two_sided_pvalue": paired_sign_flip_pvalue(difference),
                "pair_aware_wins": int((difference > 0).sum()),
                "ties": int((difference == 0).sum()),
                "pair_aware_losses": int((difference < 0).sum()),
            }
        )
    paired = pd.DataFrame(comparison_rows)

    selected_edge_rows = []
    for (replicate, method), frame in opened.groupby(["replicate", "method"], sort=True):
        valid = frame[frame.verified_BC_bracket].sort_values("proposal_rank")
        if valid.empty:
            continue
        row = valid.iloc[0].copy()
        row["edge_run_id"] = f"r{int(replicate):02d}_{method}"
        row["conditional_selection_rule"] = "highest-ranked independently verified B-C bracket"
        selected_edge_rows.append(row)
    selected_edge_pairs = pd.DataFrame(selected_edge_rows)

    OUTDIR.mkdir(parents=True)
    opened_path = OUTDIR / "frozen_proposals_with_opened_outcomes.csv"
    replicate_path = OUTDIR / "replicate_bracket_metrics.csv"
    method_path = OUTDIR / "method_summary.csv"
    paired_path = OUTDIR / "paired_method_comparisons.csv"
    edge_selection_path = OUTDIR / "selected_edge_recovery_pairs.csv"
    opened.to_csv(opened_path, index=False)
    replicate_metrics.to_csv(replicate_path, index=False)
    method_summary.to_csv(method_path, index=False)
    paired.to_csv(paired_path, index=False)
    selected_edge_pairs.to_csv(edge_selection_path, index=False)
    make_figure(OUTDIR / "matched_acquisition_bracket_yield.png", replicate_metrics, method_summary)
    summary = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "status": "frozen endpoint outcomes opened and evaluated",
        "replicate_count": int(config["initial_design"]["replicate_count"]),
        "pair_count_per_method_per_replicate": int(config["batch"]["pair_count"]),
        "endpoint_budget_per_method_per_replicate": int(config["batch"]["new_PDE_endpoint_budget"]),
        "primary_independent_unit": config["evaluation"]["primary_unit"],
        "candidate_outcomes_were_opened_only_after_proposal_hash_freeze": True,
        "frozen_proposal_manifest_sha256": sha256(proposal_manifest_path),
        "complete_label_bank_sha256": bank_manifest["complete_label_bank_sha256"],
        "method_results": method_summary.to_dict(orient="records"),
        "paired_results": paired.to_dict(orient="records"),
        "conditionally_selected_edge_recovery_pair_count": len(selected_edge_pairs),
        "edge_recovery_selection_rule": "highest-ranked valid B-C proposal for every method and matched initial-design replicate",
        "selected_edge_recovery_pairs_sha256": sha256(edge_selection_path),
        "scientific_scope": config["interpretation"],
        "next_gate": "multiclass-safe guarded refinement and prespecified exact-edge recovery",
    }
    (OUTDIR / "matched_acquisition_evaluation_summary.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8"
    )
    print(method_summary.to_string(index=False))
    print(paired.to_string(index=False))
    print(f"Matched acquisition evaluation saved in: {OUTDIR}")


if __name__ == "__main__":
    main()
