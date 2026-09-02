#!/usr/bin/env python3
"""Select one frozen multiclass pair for the next held-out round."""

import argparse
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

from multiclass_probit_gp import (
    fit_one_vs_rest,
    predict_normalized_probabilities,
    sample_joint_class_indices,
)


ROOT = Path(__file__).resolve().parent
ADAPTIVE = ROOT / "experiment_outputs" / "ml_heldout_CBA_adaptive_v1"


def sha256_file(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def verify_frozen_protocol(manifest):
    paths = {
        "protocol_hash": ROOT / "ml_heldout_CBA_protocol.json",
        "geometry_hash": ROOT / "experiment_outputs/ml_heldout_CBA_geometry_v1/ml_heldout_CBA_geometry.mat",
        "candidate_manifest_hash": ROOT / "experiment_outputs/ml_heldout_CBA_geometry_v1/heldout_candidate_manifest.csv",
        "initial_training_hash": ROOT / "experiment_outputs/ml_heldout_CBA_geometry_v1/allowed_coarse_training_labels.csv",
        "multiclass_code_hash": ROOT / "multiclass_probit_gp.py",
        "probit_code_hash": ROOT / "probit_gp.py",
        "selector_code_hash": ROOT / "select_ml_heldout_CBA_pair.py",
        "labeling_runner_hash": ROOT / "run_ml_heldout_CBA_round_labeling.m",
    }
    failed = [key for key, path in paths.items() if sha256_file(path) != manifest["hashes"][key]]
    if failed:
        raise RuntimeError("Frozen protocol source changed: " + ", ".join(failed))


def load_training(config, round_index):
    training = pd.read_csv(ROOT / config["initial_training_csv"])[
        ["lambda", "normalized_coordinate", "final_label"]
    ].copy()
    training["state_index"] = np.nan
    training["source_round"] = 0
    novelty_anchors = list(training.normalized_coordinate.astype(float))
    used_indices = set()
    for prior in range(1, round_index):
        label_file = ADAPTIVE / f"round_{prior:02d}" / "endpoint_labels.csv"
        if not label_file.is_file():
            raise FileNotFoundError(f"Round {prior} labels are required before round {round_index}")
        labels = pd.read_csv(label_file)
        if len(labels) != 2 or labels.state_index.duplicated().any():
            raise RuntimeError(f"Round {prior} endpoint-label artifact is invalid")
        used_indices |= set(labels.state_index.astype(int))
        novelty_anchors.extend(labels.normalized_coordinate.astype(float))
        known = labels[labels.final_label.isin(config["classes"])][
            ["lambda", "normalized_coordinate", "final_label", "state_index"]
        ].copy()
        known["source_round"] = prior
        training = pd.concat([training, known], ignore_index=True)
    return training, np.asarray(novelty_anchors, dtype=float), used_indices


def candidate_pairs(candidate, acquisition):
    grid = candidate.grid_index.to_numpy(int)
    first, second = np.triu_indices(len(candidate), k=1)
    separation = grid[second] - grid[first]
    keep = (separation >= int(acquisition["minimum_grid_separation"])) & (
        separation <= int(acquisition["maximum_grid_separation"])
    )
    first, second, separation = first[keep], second[keep], separation[keep]
    pair_id = [
        f"pair_{int(candidate.iloc[i].grid_index):03d}_{int(candidate.iloc[j].grid_index):03d}"
        for i, j in zip(first, second)
    ]
    return pd.DataFrame(
        {
            "pair_id": pair_id,
            "first_position": first,
            "second_position": second,
            "first_state_index": candidate.iloc[first].state_index.to_numpy(int),
            "second_state_index": candidate.iloc[second].state_index.to_numpy(int),
            "first_grid_index": candidate.iloc[first].grid_index.to_numpy(int),
            "second_grid_index": candidate.iloc[second].grid_index.to_numpy(int),
            "first_lambda": candidate.iloc[first]["lambda"].to_numpy(float),
            "second_lambda": candidate.iloc[second]["lambda"].to_numpy(float),
            "first_coordinate": candidate.iloc[first].normalized_coordinate.to_numpy(float),
            "second_coordinate": candidate.iloc[second].normalized_coordinate.to_numpy(float),
            "grid_separation": separation,
        }
    )


def save_model(path, model):
    payload = {"classes": np.asarray(model.classes)}
    for label in model.classes:
        state = model.states[label]
        prefix = f"{label}_"
        for field in ["X", "y", "lengthscales", "signal_std", "jitter", "f", "a", "sqrt_W", "L"]:
            value = getattr(state, field)
            payload[prefix + field] = np.atleast_1d(value)
    np.savez(path, **payload)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--round", type=int, required=True)
    args = parser.parse_args()
    manifest = json.loads((ADAPTIVE / "frozen_protocol_manifest.json").read_text())
    verify_frozen_protocol(manifest)
    config = json.loads((ADAPTIVE / "frozen_protocol.json").read_text())
    maximum_round = int(config["acquisition"]["round_count"])
    if args.round < 1 or args.round > maximum_round:
        raise ValueError(f"round must be between 1 and {maximum_round}")
    round_dir = ADAPTIVE / f"round_{args.round:02d}"
    if round_dir.exists():
        raise FileExistsError(f"Refusing to overwrite round selection: {round_dir}")

    candidate = pd.read_csv(ROOT / config["candidate_manifest_csv"])
    if candidate.outcome_available_to_acquisition.astype(bool).any():
        raise RuntimeError("Candidate outcome leakage detected")
    training, novelty_anchors, used_indices = load_training(config, args.round)
    available = candidate[~candidate.state_index.isin(used_indices)].copy().reset_index(drop=True)
    classes = tuple(config["classes"])
    known_training = training[training.final_label.isin(classes)].copy()
    X_train = known_training[[config["feature_column"]]].to_numpy(float)
    model = fit_one_vs_rest(
        X_train,
        known_training.final_label.to_numpy(),
        classes,
        config["optimization"],
    )
    X_candidate = available[[config["feature_column"]]].to_numpy(float)
    probabilities = predict_normalized_probabilities(model, X_candidate)
    sampled_labels, covariance_diagnostics = sample_joint_class_indices(
        model,
        X_candidate,
        config["acquisition"]["joint_sample_count"],
        config["acquisition"]["base_random_seed"] + args.round,
        config["acquisition"]["eigenvalue_floor"],
    )
    pairs = candidate_pairs(available, config["acquisition"])
    first = pairs.first_position.to_numpy(int)
    second = pairs.second_position.to_numpy(int)
    q = np.mean(sampled_labels[:, first] != sampled_labels[:, second], axis=0)
    distance = pairs.second_coordinate.to_numpy(float) - pairs.first_coordinate.to_numpy(float)
    distance_scale = float(config["acquisition"]["distance_scale_challenge_coordinates"])
    distance_penalty = np.exp(-0.5 * (distance / distance_scale) ** 2)
    midpoint = 0.5 * (
        pairs.first_coordinate.to_numpy(float) + pairs.second_coordinate.to_numpy(float)
    )
    novelty_distance = np.min(np.abs(midpoint[:, None] - novelty_anchors[None, :]), axis=1)
    novelty_scale = float(config["acquisition"]["novelty_scale_challenge_coordinates"])
    novelty_floor = float(config["acquisition"]["novelty_floor"])
    novelty_factor = novelty_floor + (1.0 - novelty_floor) * np.minimum(
        novelty_distance / novelty_scale, 1.0
    )
    score = q * distance_penalty * novelty_factor
    pairs["joint_multiclass_disagreement_probability"] = q
    pairs["distance_penalty"] = distance_penalty
    pairs["midpoint_coordinate"] = midpoint
    pairs["midpoint_novelty_distance"] = novelty_distance
    pairs["novelty_factor"] = novelty_factor
    pairs["acquisition_score"] = score
    for class_index, label in enumerate(classes):
        pairs[f"first_probability_{label}"] = probabilities[first, class_index]
        pairs[f"second_probability_{label}"] = probabilities[second, class_index]
    pairs = pairs.sort_values(
        ["acquisition_score", "joint_multiclass_disagreement_probability", "grid_separation", "pair_id"],
        ascending=[False, False, True, True],
    ).reset_index(drop=True)
    pairs.insert(1, "score_rank", np.arange(1, len(pairs) + 1))
    selection = pairs.iloc[[0]].copy()
    selection.insert(0, "round", args.round)

    round_dir.mkdir(parents=True)
    pairs.to_csv(round_dir / "pair_ranking.csv", index=False)
    selection.to_csv(round_dir / "selected_pair.csv", index=False)
    known_training.to_csv(round_dir / "training_snapshot.csv", index=False)
    save_model(round_dir / "multiclass_probit_gp_model.npz", model)
    diagnostics = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "round": args.round,
        "training_count": int(len(known_training)),
        "training_class_counts": known_training.final_label.value_counts().to_dict(),
        "available_candidate_count": int(len(available)),
        "candidate_pair_count": int(len(pairs)),
        "selected_pair_id": str(selection.pair_id.iloc[0]),
        "selected_pair_score": float(selection.acquisition_score.iloc[0]),
        "selected_pair_q": float(selection.joint_multiclass_disagreement_probability.iloc[0]),
        "covariance_diagnostics": covariance_diagnostics,
        "hyperparameters": {
            label: {
                "lengthscale": float(model.states[label].lengthscales[0]),
                "signal_std": float(model.states[label].signal_std),
            }
            for label in classes
        },
        "new_heldout_outcomes_used_for_this_selection": args.round > 1,
        "outcome_use_rule": "only completed earlier adaptive rounds; no pilot-refinement outcome file",
    }
    (round_dir / "selection_diagnostics.json").write_text(json.dumps(diagnostics, indent=2))
    selection_manifest = {
        "round": args.round,
        "selection_frozen_before_round_labels": True,
        "selected_pair_hash": sha256_file(round_dir / "selected_pair.csv"),
        "training_snapshot_hash": sha256_file(round_dir / "training_snapshot.csv"),
        "model_hash": sha256_file(round_dir / "multiclass_probit_gp_model.npz"),
        "selected_state_indices": [
            int(selection.first_state_index.iloc[0]),
            int(selection.second_state_index.iloc[0]),
        ],
        "selected_lambdas": [
            float(selection.first_lambda.iloc[0]),
            float(selection.second_lambda.iloc[0]),
        ],
        "labels_opened": False,
    }
    (round_dir / "selection_manifest.json").write_text(json.dumps(selection_manifest, indent=2))

    figure, axis = plt.subplots(figsize=(8, 4.2), constrained_layout=True)
    for class_index, label in enumerate(classes):
        axis.plot(available["lambda"], probabilities[:, class_index], label=f"p({label})")
    for value in selection[["first_lambda", "second_lambda"]].iloc[0]:
        axis.axvline(value, color="black", linestyle="--")
    axis.scatter(known_training["lambda"], np.full(len(known_training), -0.04), c="0.3", s=14)
    axis.set_ylim(-0.08, 1.02)
    axis.set_xlabel("lambda")
    axis.set_ylabel("normalized one-vs-rest probability")
    axis.set_title(f"Frozen held-out pair selection, round {args.round}")
    axis.legend(ncol=3)
    figure.savefig(round_dir / "selection_probabilities.png", dpi=300)
    plt.close(figure)
    print(json.dumps(diagnostics, indent=2))
    print(selection.to_string(index=False))


if __name__ == "__main__":
    main()
