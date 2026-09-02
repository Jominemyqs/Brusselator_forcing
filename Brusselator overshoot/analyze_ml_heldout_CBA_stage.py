#!/usr/bin/env python3
"""Consolidate the frozen held-out test and post-budget guarded refinement."""

import json
import os
from datetime import datetime, timezone
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/private/tmp/brusselator_matplotlib")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parent
EVALUATION = ROOT / "experiment_outputs" / "ml_heldout_CBA_evaluation_v1"
REFINEMENT = ROOT / "experiment_outputs" / "ml_heldout_CBA_postbudget_refinement_v1"
GEOMETRY = ROOT / "experiment_outputs" / "ml_heldout_CBA_geometry_v1"
OUTDIR = ROOT / "experiment_outputs" / "ml_heldout_CBA_stage_summary_v1"


def main():
    if OUTDIR.exists():
        raise FileExistsError(f"Refusing to overwrite held-out stage summary: {OUTDIR}")
    fixed = json.loads((EVALUATION / "heldout_evaluation.json").read_text())
    adaptive = pd.read_csv(EVALUATION / "all_adaptive_labels.csv")
    midpoint = pd.read_csv(REFINEMENT / "postbudget_midpoint_cases.csv")
    brackets = pd.read_csv(REFINEMENT / "discovered_adjacent_brackets.csv")
    coarse = pd.read_csv(GEOMETRY / "allowed_coarse_training_labels.csv")
    if fixed["primary_B_detection_success"] or fixed["new_PDE_runs"] != 8:
        raise RuntimeError("The frozen-budget result changed")
    if set(brackets.bracket_id) != {"CB", "BA"}:
        raise RuntimeError("Post-budget refinement did not produce both expected brackets")
    postbudget_B = midpoint[midpoint.outcome == "B"]
    if len(postbudget_B) != 1:
        raise RuntimeError("Expected exactly one post-budget B discovery state")

    fixed_runtime = float(adaptive.runtime_seconds.sum())
    refinement_runtime = float(midpoint.runtime_seconds.sum())
    summary = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "frozen_eight_run_test": {
            "primary_B_detection_success": False,
            "labels": adaptive.final_label.value_counts().to_dict(),
            "different_label_pairs": int(fixed["different_label_selected_pairs"]),
            "final_CA_interval": fixed["final_selected_pair"],
            "recorded_runtime_seconds": fixed_runtime,
        },
        "postbudget_guarded_refinement": {
            "additional_PDE_runs": int(len(midpoint)),
            "outcome_sequence": midpoint.outcome.tolist(),
            "B_lambda": float(postbudget_B["lambda"].iloc[0]),
            "recorded_runtime_seconds": refinement_runtime,
            "brackets": brackets[
                ["bracket_id", "left_lambda", "left_outcome", "right_lambda", "right_outcome", "width"]
            ].to_dict(orient="records"),
        },
        "total_new_PDE_runs_through_B_discovery": int(len(adaptive) + len(midpoint)),
        "total_recorded_runtime_seconds": fixed_runtime + refinement_runtime,
        "scientific_interpretation": (
            "the frozen ML policy did not directly sample the narrow B basin under its fixed budget; "
            "it localized a C/A interval containing B, after which four multiclass-guarded midpoint "
            "runs recovered the C-B-A structure and produced two edge-tracking brackets"
        ),
        "next_dynamical_tests": [
            "edge-track the C-B bracket and test recovery of the independently known E_BC orbit",
            "edge-track the B-A bracket without assuming its organizing invariant object",
        ],
    }

    OUTDIR.mkdir(parents=True)
    with (OUTDIR / "heldout_stage_summary.json").open("w") as handle:
        json.dump(summary, handle, indent=2)
    brackets.to_csv(OUTDIR / "edge_tracking_targets.csv", index=False)

    colors = {"A": "tab:blue", "B": "tab:orange", "C": "gold"}
    figure, axes = plt.subplots(1, 3, figsize=(14, 4.4), constrained_layout=True)
    local_coarse = coarse[(coarse["lambda"] >= 0.275) & (coarse["lambda"] <= 0.3)]
    ymap = {"C": 1, "B": 2, "A": 3}
    for label in ["C", "B", "A"]:
        group = adaptive[adaptive.final_label == label]
        axes[0].scatter(group["lambda"], [ymap[label]] * len(group), color=colors[label], s=55, label=f"adaptive {label}")
    axes[0].scatter(local_coarse["lambda"], local_coarse.final_label.map(ymap), marker="s", color="0.3", label="coarse endpoints")
    axes[0].set_yticks([1, 2, 3], ["C", "B", "A"])
    axes[0].set_xlabel("lambda")
    axes[0].set_title("Frozen eight-run outcomes")
    axes[0].legend(fontsize=8)

    for label in ["C", "B", "A"]:
        group = midpoint[midpoint.outcome == label]
        axes[1].scatter(group["lambda"], [ymap[label]] * len(group), color=colors[label], s=75, label=label)
    for row in brackets.itertuples():
        axes[1].axvspan(row.left_lambda, row.right_lambda, alpha=0.18, label=f"{row.bracket_id} bracket")
    axes[1].set_yticks([1, 2, 3], ["C", "B", "A"])
    axes[1].set_xlabel("lambda")
    axes[1].set_title("Post-budget guarded refinement")
    axes[1].legend(fontsize=8)

    decision = list(adaptive.decision_time) + list(midpoint.decision_time)
    bar_colors = [colors.get(label, "0.4") for label in list(adaptive.final_label) + list(midpoint.outcome)]
    axes[2].bar(np.arange(1, len(decision) + 1), decision, color=bar_colors)
    axes[2].axvline(8.5, color="black", linestyle="--")
    axes[2].set_xlabel("new PDE run")
    axes[2].set_ylabel("decision time")
    axes[2].set_title("Fixed budget | post-budget")
    figure.suptitle("Held-out C--B--A challenge: localization followed by guarded discovery")
    figure.savefig(OUTDIR / "heldout_CBA_stage.png", dpi=300)
    plt.close(figure)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
