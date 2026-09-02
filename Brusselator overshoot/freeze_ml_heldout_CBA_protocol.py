#!/usr/bin/env python3
"""Freeze the multiclass held-out acquisition policy before new PDE labels."""

import hashlib
import json
import shutil
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parent
CONFIG = ROOT / "ml_heldout_CBA_protocol.json"


def sha256_file(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main():
    config = json.loads(CONFIG.read_text(encoding="utf-8"))
    outdir = ROOT / config["output_directory"]
    if outdir.exists():
        raise FileExistsError(f"Refusing to overwrite frozen protocol: {outdir}")
    candidates = pd.read_csv(ROOT / config["candidate_manifest_csv"])
    training = pd.read_csv(ROOT / config["initial_training_csv"])
    if candidates.outcome_available_to_acquisition.astype(bool).any():
        raise RuntimeError("A held-out candidate outcome is exposed")
    if set(training.outcome_source) != {"coarse_scan_only"}:
        raise RuntimeError("Training data are not restricted to the coarse scan")
    if set(training.final_label) != {"A", "B", "C"}:
        raise RuntimeError("Initial multiclass training taxonomy is incomplete")

    tracked = {
        "protocol_hash": CONFIG,
        "geometry_hash": ROOT / config["geometry_mat"],
        "candidate_manifest_hash": ROOT / config["candidate_manifest_csv"],
        "initial_training_hash": ROOT / config["initial_training_csv"],
        "multiclass_code_hash": ROOT / "multiclass_probit_gp.py",
        "probit_code_hash": ROOT / "probit_gp.py",
        "selector_code_hash": ROOT / "select_ml_heldout_CBA_pair.py",
        "labeling_runner_hash": ROOT / "run_ml_heldout_CBA_round_labeling.m",
    }
    missing = [str(path) for path in tracked.values() if not path.is_file()]
    if missing:
        raise FileNotFoundError("Protocol source is missing: " + ", ".join(missing))
    manifest = {
        "schema_version": "1.0",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "protocol_frozen_before_new_labels": True,
        "new_heldout_outcomes_used": False,
        "pairwise_refinement_outcomes_used": False,
        "candidate_count": int(len(candidates)),
        "initial_training_count": int(len(training)),
        "round_count": int(config["acquisition"]["round_count"]),
        "total_new_PDE_budget": int(config["acquisition"]["total_new_PDE_budget"]),
        "hashes": {key: sha256_file(path) for key, path in tracked.items()},
        "status": "frozen and unopened",
    }
    outdir.mkdir(parents=True)
    shutil.copy2(CONFIG, outdir / "frozen_protocol.json")
    (outdir / "frozen_protocol_manifest.json").write_text(
        json.dumps(manifest, indent=2), encoding="utf-8"
    )
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
