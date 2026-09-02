#!/usr/bin/env python3
"""Freeze nested matched designs and implementation hashes for replay robustness."""

from __future__ import annotations

import hashlib
import json
import shutil
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parent
CONFIG_FILE = ROOT / "ml_bc_acquisition_robustness_config.json"
IMPLEMENTATION_FILES = (
    "freeze_ml_BC_acquisition_robustness_protocol.py",
    "run_ml_BC_acquisition_robustness_proposals.py",
    "evaluate_ml_BC_acquisition_robustness.py",
    "probit_gp_kernels.py",
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


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


def main() -> None:
    config = json.loads(CONFIG_FILE.read_text(encoding="utf-8"))
    outdir = ROOT / config["protocol_output_directory"]
    if outdir.exists():
        raise FileExistsError(f"Refusing to overwrite frozen robustness protocol: {outdir}")
    for relative in IMPLEMENTATION_FILES:
        if not (ROOT / relative).exists():
            raise FileNotFoundError(f"Freeze requires completed implementation file: {relative}")

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

    bank_path = ROOT / config["complete_label_bank_csv"]
    bank_manifest_path = ROOT / config["complete_label_bank_manifest_json"]
    bank_manifest = json.loads(bank_manifest_path.read_text(encoding="utf-8"))
    if sha256(bank_path) != bank_manifest["complete_label_bank_sha256"]:
        raise ValueError("Complete label bank hash changed")
    bank = pd.read_csv(bank_path, usecols=["state_index", "final_label"])
    label_lookup = dict(zip(bank.state_index.astype(int), bank.final_label.astype(str)))

    design_cfg = config["matched_designs"]
    anchors = [int(item["state_index"]) for item in design_cfg["fixed_anchors"]]
    for item in design_cfg["fixed_anchors"]:
        state_index = int(item["state_index"])
        if label_lookup.get(state_index) != str(item["required_label"]):
            raise ValueError(f"Fixed anchor {state_index} has the wrong saved label")

    coordinates = geometry[config["model"]["feature_columns"]].to_numpy(float)
    state_indices = geometry.state_index.to_numpy(int)
    maximum_count = int(design_cfg["maximum_training_count"])
    prefixes = [int(value) for value in design_cfg["training_count_prefixes"]]
    rows = []
    designs_by_prefix = {prefix: set() for prefix in prefixes}
    for replicate in range(1, int(design_cfg["replicate_count"]) + 1):
        seed = int(design_cfg["base_random_seed"]) + replicate - 1
        selected = randomized_maximin_design(
            coordinates,
            state_indices,
            anchors,
            maximum_count,
            int(design_cfg["randomized_top_k"]),
            seed,
        )
        for prefix in prefixes:
            key = tuple(sorted(selected[:prefix]))
            if key in designs_by_prefix[prefix]:
                raise ValueError(f"Duplicate matched design at prefix n0={prefix}")
            designs_by_prefix[prefix].add(key)
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
                    "in_n15_prefix": position <= 15,
                    "in_n25_prefix": position <= 25,
                    "in_n35_prefix": position <= 35,
                }
            )
    designs = pd.DataFrame(rows)
    if not np.all(designs.groupby("replicate").size() == maximum_count):
        raise ValueError("A matched design has the wrong maximum size")

    outdir.mkdir(parents=True)
    shutil.copy2(CONFIG_FILE, outdir / CONFIG_FILE.name)
    designs_path = outdir / "nested_matched_initial_designs.csv"
    designs.to_csv(designs_path, index=False)
    scenario_path = outdir / "frozen_scenarios.csv"
    pd.DataFrame(config["scenarios"]).to_csv(scenario_path, index=False)

    implementation_hashes = {relative: sha256(ROOT / relative) for relative in IMPLEMENTATION_FILES}
    manifest = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "protocol_frozen": True,
        "protocol_version": config["experiment_name"],
        "configuration_sha256": sha256(CONFIG_FILE),
        "geometry_manifest_sha256": sha256(geometry_path),
        "complete_label_bank_sha256": bank_manifest["complete_label_bank_sha256"],
        "complete_label_bank_manifest_sha256": sha256(bank_manifest_path),
        "nested_matched_initial_designs_sha256": sha256(designs_path),
        "frozen_scenarios_sha256": sha256(scenario_path),
        "implementation_sha256": implementation_hashes,
        "replicate_count": int(design_cfg["replicate_count"]),
        "maximum_training_count": maximum_count,
        "training_count_prefixes": prefixes,
        "nonanchor_outcomes_used_for_design": False,
        "fixed_anchor_labels_used": True,
        "new_PDE_integrations": 0,
        "masked_replay": True,
        "probability_index_audit_correction": config["audit_correction"],
        "scientific_scope": config["evaluation_scope"]["claim"],
        "excluded_claims": config["evaluation_scope"]["excluded_claims"],
    }
    manifest_path = outdir / "frozen_robustness_protocol_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(json.dumps(manifest, indent=2))
    print(f"Frozen robustness protocol saved in: {outdir}")


if __name__ == "__main__":
    main()
