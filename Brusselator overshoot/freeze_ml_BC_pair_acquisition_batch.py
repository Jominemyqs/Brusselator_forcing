#!/usr/bin/env python3
"""Freeze the robust nominal four-pair development batch before PDE labels."""

import hashlib
import json
import shutil
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parent
SOURCE = ROOT / "experiment_outputs" / "ml_BC_pair_acquisition_sensitivity_v1"
OUTDIR = ROOT / "experiment_outputs" / "ml_BC_pair_acquisition_frozen_v1"


def sha256_file(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main():
    if OUTDIR.exists():
        raise FileExistsError(f"Refusing to overwrite frozen pair batch: {OUTDIR}")
    summary = json.loads((SOURCE / "sensitivity_summary.json").read_text(encoding="utf-8"))
    if summary["new_PDE_outcomes_used"] or summary["primary_batch_frozen"]:
        raise RuntimeError("Sensitivity source is not in the required outcome-blind state")
    if summary["minimum_top_k_jaccard"] < 0.8 or summary["minimum_four_pair_batch_intersection"] < 4:
        raise RuntimeError("Training-only ranking robustness is insufficient to freeze the batch")

    batch = pd.read_csv(SOURCE / "nominal_saved_four_pair_batch.csv")
    if len(batch) != 4:
        raise ValueError("The development batch must contain four pairs")
    endpoints = []
    for row in batch.itertuples(index=False):
        for side in ["first", "second"]:
            endpoints.append(
                {
                    "state_index": int(getattr(row, f"{side}_state_index")),
                    "sample_id": getattr(row, f"{side}_sample_id"),
                    "normalized_alpha_1": float(getattr(row, f"{side}_normalized_alpha_1")),
                    "normalized_alpha_2": float(getattr(row, f"{side}_normalized_alpha_2")),
                    "pair_batch_step": int(row.batch_step),
                    "pair_id": row.pair_id,
                    "pair_side": side,
                }
            )
    endpoint_frame = pd.DataFrame(endpoints).sort_values("state_index")
    if len(endpoint_frame) != 8 or endpoint_frame.state_index.nunique() != 8:
        raise ValueError("The four frozen pairs must use eight distinct endpoints")

    OUTDIR.mkdir(parents=True)
    batch.to_csv(OUTDIR / "frozen_pair_batch.csv", index=False)
    endpoint_frame.to_csv(OUTDIR / "frozen_endpoint_manifest.csv", index=False)
    shutil.copy2(SOURCE / "ml_bc_pair_acquisition_sensitivity_config.json", OUTDIR)
    shutil.copy2(SOURCE / "pair_ranking_sensitivity.csv", OUTDIR)
    manifest = {
        "schema_version": "1.0",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "experiment_name": "ml_BC_pair_acquisition_frozen_v1",
        "pair_batch_frozen": True,
        "endpoint_outcomes_used": False,
        "fixed_test_outcomes_used_for_selection": False,
        "revealed_outcomes_used_for_selection": False,
        "revealed_state_indices_used_for_exclusion_only": True,
        "primary_fit": "nominal_saved",
        "pair_count": 4,
        "endpoint_count": 8,
        "new_PDE_endpoint_budget": 8,
        "minimum_top_20_jaccard_across_sensitivity_fits": summary["minimum_top_k_jaccard"],
        "four_pair_batch_intersection_across_every_fit": summary["minimum_four_pair_batch_intersection"],
        "frozen_pair_batch_hash": sha256_file(OUTDIR / "frozen_pair_batch.csv"),
        "frozen_endpoint_manifest_hash": sha256_file(OUTDIR / "frozen_endpoint_manifest.csv"),
        "nominal_model_hash": sha256_file(
            ROOT / "experiment_outputs" / "ml_BC_probabilistic_boundary_v1" / "probit_gp_model.npz"
        ),
        "sensitivity_summary_hash": sha256_file(SOURCE / "sensitivity_summary.json"),
        "exact_edge_orbit_used": False,
        "held_out_CA_region_used": False,
        "phase": "prospective endpoint labeling inside a retrospective B-C method-development benchmark",
    }
    with (OUTDIR / "frozen_pair_manifest.json").open("w", encoding="utf-8") as handle:
        json.dump(manifest, handle, indent=2)
    print(batch[[
        "batch_step", "pair_id", "first_state_index", "second_state_index",
        "joint_opposite_sign_probability", "acquisition_score"
    ]].to_string(index=False))
    print(json.dumps(manifest, indent=2))
    print(f"Frozen pair batch saved in: {OUTDIR}")


if __name__ == "__main__":
    main()
