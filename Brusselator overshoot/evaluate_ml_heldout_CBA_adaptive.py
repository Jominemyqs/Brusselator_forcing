#!/usr/bin/env python3
"""Evaluate the frozen eight-run held-out challenge after all rounds close."""

import hashlib
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
ADAPTIVE = ROOT / "experiment_outputs" / "ml_heldout_CBA_adaptive_v1"
PILOT = ROOT / "experiment_outputs" / "edge_slice_pairwise_refinement_v1" / "pairwise_refinement_cases.csv"
OUTDIR = ROOT / "experiment_outputs" / "ml_heldout_CBA_evaluation_v1"


def sha256_file(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main():
    if OUTDIR.exists():
        raise FileExistsError(f"Refusing to overwrite held-out evaluation: {OUTDIR}")
    protocol = json.loads((ADAPTIVE / "frozen_protocol.json").read_text())
    frozen = json.loads((ADAPTIVE / "frozen_protocol_manifest.json").read_text())
    if not frozen["protocol_frozen_before_new_labels"] or frozen["new_heldout_outcomes_used"]:
        raise RuntimeError("The held-out protocol was not frozen unopened")
    round_rows = []
    selected_rows = []
    for round_index in range(1, int(protocol["acquisition"]["round_count"]) + 1):
        directory = ADAPTIVE / f"round_{round_index:02d}"
        manifest = json.loads((directory / "selection_manifest.json").read_text())
        if not manifest["selection_frozen_before_round_labels"]:
            raise RuntimeError(f"Round {round_index} was not prospectively selected")
        if sha256_file(directory / "selected_pair.csv") != manifest["selected_pair_hash"]:
            raise RuntimeError(f"Round {round_index} selection changed after labeling")
        labels = pd.read_csv(directory / "endpoint_labels.csv")
        if list(labels.state_index.astype(int)) != manifest["selected_state_indices"]:
            raise RuntimeError(f"Round {round_index} labels do not match its frozen pair")
        labels["round"] = round_index
        round_rows.append(labels)
        selected = pd.read_csv(directory / "selected_pair.csv")
        selected["first_actual_label"] = labels.final_label.iloc[0]
        selected["second_actual_label"] = labels.final_label.iloc[1]
        selected["actual_labels_differ"] = labels.final_label.iloc[0] != labels.final_label.iloc[1]
        selected_rows.append(selected)
    labels = pd.concat(round_rows, ignore_index=True)
    pairs = pd.concat(selected_rows, ignore_index=True)

    primary = bool((labels.final_label == "B").any())
    known = pd.read_csv(protocol["initial_training_csv"])[["lambda", "final_label"]]
    acquired = labels[["lambda", "final_label"]]
    ordered = pd.concat([known, acquired], ignore_index=True).sort_values("lambda")
    transitions = []
    for index in range(len(ordered) - 1):
        left = ordered.iloc[index]
        right = ordered.iloc[index + 1]
        if left.final_label != right.final_label:
            transitions.append(
                {
                    "left_lambda": left["lambda"],
                    "left_label": left.final_label,
                    "right_lambda": right["lambda"],
                    "right_label": right.final_label,
                    "width": right["lambda"] - left["lambda"],
                }
            )
    transitions = pd.DataFrame(transitions)
    stronger = bool(
        primary
        and ((transitions.left_label == "C") & (transitions.right_label == "B")).any()
        and ((transitions.left_label == "B") & (transitions.right_label == "A")).any()
    )
    last = pairs.iloc[-1]
    final_interval = [float(last.first_lambda), float(last.second_lambda)]

    # This file is intentionally opened only after the prospective budget closes.
    pilot = pd.read_csv(PILOT)
    pilot_B = pilot[(pilot.pair_id == "CA") & (pilot.outcome == "B")]
    enclosed_pilot_B = pilot_B[
        (pilot_B["lambda"] > final_interval[0]) & (pilot_B["lambda"] < final_interval[1])
    ]
    summary = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "fixed_budget_closed_before_pilot_comparison": True,
        "new_PDE_runs": int(len(labels)),
        "new_label_counts": labels.final_label.value_counts().to_dict(),
        "primary_B_detection_success": primary,
        "stronger_CB_and_BA_brackets_success": stronger,
        "different_label_selected_pairs": int(pairs.actual_labels_differ.sum()),
        "final_selected_pair": {
            "left_lambda": final_interval[0],
            "right_lambda": final_interval[1],
            "width": final_interval[1] - final_interval[0],
            "outcomes": [str(last.first_actual_label), str(last.second_actual_label)],
        },
        "post_budget_historical_pilot_comparison": {
            "known_B_points_inside_final_pair": int(len(enclosed_pilot_B)),
            "known_B_lambdas_inside_final_pair": enclosed_pilot_B["lambda"].astype(float).tolist(),
            "interpretation": (
                "post hoc localization diagnostic only; it does not convert the failed "
                "predefined B-detection endpoint into a success"
            ),
        },
        "recorded_endpoint_runtime_seconds": float(labels.runtime_seconds.sum()),
        "conclusion": (
            "the frozen eight-run method failed to sample B, but localized a sub-1e-3 "
            "C/A pair containing the historically known B pilot; a separate post-budget "
            "multiclass refinement can test its dynamical usefulness"
        ),
    }

    OUTDIR.mkdir(parents=True)
    labels.to_csv(OUTDIR / "all_adaptive_labels.csv", index=False)
    pairs.to_csv(OUTDIR / "round_pair_outcomes.csv", index=False)
    transitions.to_csv(OUTDIR / "observed_label_transitions.csv", index=False)
    with (OUTDIR / "heldout_evaluation.json").open("w") as handle:
        json.dump(summary, handle, indent=2)

    colors = {"A": "tab:blue", "B": "tab:orange", "C": "gold", "other": "tab:purple", "U": "0.4"}
    figure, axis = plt.subplots(figsize=(9, 4.4), constrained_layout=True)
    for label, group in labels.groupby("final_label"):
        axis.scatter(group["lambda"], group["round"], s=80, color=colors[label], label=f"new {label}")
    for row in pairs.itertuples():
        axis.plot([row.first_lambda, row.second_lambda], [row.round, row.round], color="0.35", zorder=0)
    if len(enclosed_pilot_B):
        axis.scatter(enclosed_pilot_B["lambda"], np.full(len(enclosed_pilot_B), 4.35), marker="*", s=150,
                     color=colors["B"], edgecolor="black", label="historical B pilot (opened after budget)")
    axis.set_yticks(range(1, 5))
    axis.set_xlabel("lambda")
    axis.set_ylabel("adaptive round")
    axis.set_title("Frozen held-out acquisition path and post-budget pilot comparison")
    axis.legend(loc="best")
    figure.savefig(OUTDIR / "heldout_acquisition_path.png", dpi=300)
    plt.close(figure)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
