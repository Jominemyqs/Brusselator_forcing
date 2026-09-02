#!/usr/bin/env python3
"""Freeze matched initial designs and the remaining-label manifest.

The design construction reads geometry and the two explicitly declared anchor
labels. No other outcome is used. The resulting protocol must exist before the
remaining frozen-PDE labels are generated.
"""

from __future__ import annotations

import hashlib
import json
import shutil
from datetime import datetime, timezone
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parent
CONFIG_FILE = ROOT / "ml_bc_matched_benchmark_protocol_config.json"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_existing_indices(config: dict) -> tuple[set[int], dict[int, str], list[dict]]:
    indices: set[int] = set()
    labels: dict[int, str] = {}
    sources = []
    for relative in config["existing_label_sources"]:
        path = ROOT / relative
        frame = pd.read_csv(path, usecols=["state_index", "final_label"])
        for row in frame.itertuples(index=False):
            state_index = int(row.state_index)
            label = str(row.final_label)
            if state_index in labels and labels[state_index] != label:
                raise ValueError(f"Conflicting saved labels for state {state_index}")
            labels[state_index] = label
            indices.add(state_index)
        sources.append(
            {
                "path": relative,
                "sha256": sha256(path),
                "row_count": int(len(frame)),
            }
        )
    return indices, labels, sources


def randomized_maximin_design(
    coordinates: np.ndarray,
    state_indices: np.ndarray,
    anchors: list[int],
    count: int,
    top_k: int,
    seed: int,
) -> list[int]:
    position = {int(state): pos for pos, state in enumerate(state_indices)}
    selected = list(anchors)
    available = set(map(int, state_indices)) - set(selected)
    rng = np.random.default_rng(seed)
    while len(selected) < count:
        candidates = np.array(sorted(available), dtype=int)
        candidate_pos = np.array([position[value] for value in candidates])
        selected_pos = np.array([position[value] for value in selected])
        delta = coordinates[candidate_pos, None, :] - coordinates[selected_pos, :][None, :, :]
        minimum_distance = np.sqrt(np.sum(delta * delta, axis=2)).min(axis=1)
        ordering = np.lexsort((candidates, -minimum_distance))
        shortlist = candidates[ordering[: min(top_k, len(ordering))]]
        chosen = int(shortlist[rng.integers(0, len(shortlist))])
        selected.append(chosen)
        available.remove(chosen)
    return selected


def make_design_figure(path: Path, geometry: pd.DataFrame, designs: pd.DataFrame) -> None:
    figure, axes = plt.subplots(3, 4, figsize=(12, 9), constrained_layout=True)
    for replicate, axis in enumerate(axes.flat, start=1):
        chosen = designs[designs.replicate == replicate]
        axis.scatter(
            geometry.normalized_alpha_1,
            geometry.normalized_alpha_2,
            s=8,
            color="0.85",
        )
        nonanchor = chosen[~chosen.is_fixed_anchor]
        anchor = chosen[chosen.is_fixed_anchor]
        axis.scatter(
            nonanchor.normalized_alpha_1,
            nonanchor.normalized_alpha_2,
            s=24,
            color="#0072B2",
        )
        axis.scatter(
            anchor.normalized_alpha_1,
            anchor.normalized_alpha_2,
            s=52,
            marker="*",
            color="#D55E00",
        )
        axis.set_title(f"design {replicate:02d}")
        axis.set_aspect("equal")
        axis.set_xlim(-1.05, 1.05)
        axis.set_ylim(-1.05, 1.05)
    figure.suptitle("Frozen outcome-blind matched initial designs")
    figure.supxlabel(r"normalized $\alpha_1$")
    figure.supylabel(r"normalized $\alpha_2$")
    figure.savefig(path, dpi=250)
    plt.close(figure)


def main() -> None:
    config = json.loads(CONFIG_FILE.read_text(encoding="utf-8"))
    outdir = ROOT / config["output_root"] / config["experiment_name"]
    if outdir.exists():
        raise FileExistsError(f"Refusing to overwrite frozen protocol: {outdir}")

    geometry_path = ROOT / config["geometry_manifest_csv"]
    geometry = pd.read_csv(geometry_path).sort_values("state_index").reset_index(drop=True)
    gate = config["geometry_gate"]
    if len(geometry) != int(gate["expected_state_count"]):
        raise ValueError("Frozen plane state count changed")
    if geometry.state_index.tolist() != list(range(1, len(geometry) + 1)):
        raise ValueError("Frozen plane state indices are no longer contiguous")
    if geometry.row_index.nunique() != int(gate["expected_grid_rows"]) or geometry.column_index.nunique() != int(
        gate["expected_grid_columns"]
    ):
        raise ValueError("Frozen plane dimensions changed")

    known_indices, known_labels, label_sources = read_existing_indices(config)
    design_cfg = config["initial_design"]
    anchors = [int(item["state_index"]) for item in design_cfg["fixed_anchors"]]
    for item in design_cfg["fixed_anchors"]:
        state_index = int(item["state_index"])
        required = str(item["required_label"])
        if known_labels.get(state_index) != required:
            raise ValueError(f"Required anchor {state_index} is not saved as {required}")

    coordinates = geometry[["normalized_alpha_1", "normalized_alpha_2"]].to_numpy(float)
    state_indices = geometry.state_index.to_numpy(int)
    rows = []
    seen_designs: set[tuple[int, ...]] = set()
    for replicate in range(1, int(design_cfg["replicate_count"]) + 1):
        seed = int(design_cfg["base_random_seed"]) + replicate - 1
        selected = randomized_maximin_design(
            coordinates,
            state_indices,
            anchors,
            int(design_cfg["training_count"]),
            int(design_cfg["randomized_top_k"]),
            seed,
        )
        key = tuple(sorted(selected))
        if key in seen_designs:
            raise ValueError("Randomized maximin produced duplicate initial designs")
        seen_designs.add(key)
        for position, state_index in enumerate(selected, start=1):
            item = geometry.loc[geometry.state_index == state_index].iloc[0]
            rows.append(
                {
                    "replicate": replicate,
                    "design_seed": seed,
                    "design_position": position,
                    "state_index": state_index,
                    "sample_id": item.sample_id,
                    "normalized_alpha_1": item.normalized_alpha_1,
                    "normalized_alpha_2": item.normalized_alpha_2,
                    "is_fixed_anchor": state_index in anchors,
                }
            )
    designs = pd.DataFrame(rows)
    if not np.all(designs.groupby("replicate").size() == int(design_cfg["training_count"])):
        raise ValueError("A frozen initial design has the wrong size")

    missing = geometry[~geometry.state_index.isin(known_indices)].copy()
    missing.insert(0, "manifest_position", np.arange(1, len(missing) + 1))
    known = geometry[geometry.state_index.isin(known_indices)].copy()
    known["saved_label"] = known.state_index.map(known_labels)

    outdir.mkdir(parents=True)
    shutil.copy2(CONFIG_FILE, outdir / CONFIG_FILE.name)
    designs_path = outdir / "matched_initial_designs.csv"
    missing_path = outdir / "missing_label_manifest.csv"
    known_path = outdir / "existing_label_index.csv"
    designs.to_csv(designs_path, index=False)
    missing.to_csv(missing_path, index=False)
    known.to_csv(known_path, index=False)
    make_design_figure(outdir / "matched_initial_designs.png", geometry, designs)

    manifest = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "protocol_frozen": True,
        "protocol_version": config["experiment_name"],
        "geometry_id": gate["expected_geometry_id"],
        "geometry_manifest_sha256": sha256(geometry_path),
        "configuration_sha256": sha256(CONFIG_FILE),
        "matched_initial_designs_sha256": sha256(designs_path),
        "missing_label_manifest_sha256": sha256(missing_path),
        "replicate_count": int(design_cfg["replicate_count"]),
        "training_count_per_replicate": int(design_cfg["training_count"]),
        "pair_count_per_method": int(config["batch"]["pair_count"]),
        "endpoint_budget_per_method": int(config["batch"]["new_PDE_endpoint_budget"]),
        "existing_labeled_state_count": len(known_indices),
        "missing_label_state_count": len(missing),
        "labels_used_for_design": False,
        "fixed_anchor_labels_used": True,
        "nonanchor_labels_used": False,
        "candidate_outcomes_used": False,
        "exact_edge_orbit_used": False,
        "heldout_CBA_region_used": False,
        "existing_label_sources": label_sources,
        "interpretation": config["interpretation"],
    }
    manifest_path = outdir / "frozen_protocol_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(json.dumps(manifest, indent=2))
    print(f"Frozen matched benchmark protocol saved in: {outdir}")


if __name__ == "__main__":
    main()
