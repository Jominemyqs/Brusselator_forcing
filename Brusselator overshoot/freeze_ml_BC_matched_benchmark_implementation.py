#!/usr/bin/env python3
"""Hash-freeze the matched benchmark implementation before label-bank completion."""

from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parent
OUTFILE = (
    ROOT
    / "experiment_outputs"
    / "ml_BC_matched_benchmark_protocol_v1"
    / "frozen_implementation_manifest_v7.json"
)
FILES = [
    "ml_bc_matched_benchmark_protocol_config.json",
    "freeze_ml_BC_matched_benchmark_protocol.py",
    "brusselator_ml_matched_label_bank_config.m",
    "run_ml_BC_matched_label_bank.m",
    "run_ml_BC_subset_labeling.m",
    "assemble_ml_BC_matched_label_bank.py",
    "freeze_ml_BC_matched_acquisition_proposals.py",
    "evaluate_ml_BC_matched_acquisition.py",
    "prepare_ml_BC_matched_unique_edge_pairs.py",
    "brusselator_ml_matched_edge_recovery_config.m",
    "run_ml_BC_matched_edge_recovery.m",
    "analyze_ml_BC_matched_edge_recovery.py",
    "run_ml_BC_matched_benchmark_pipeline.m",
    "freeze_ml_BC_matched_label_bank_shards.py",
    "brusselator_ml_matched_label_bank_shard_config.m",
    "run_ml_BC_matched_label_bank_shard.m",
    "run_ml_BC_matched_sharded_pipeline.sh",
]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    if OUTFILE.exists():
        raise FileExistsError(f"Refusing to overwrite frozen implementation: {OUTFILE}")
    payload = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "implementation_frozen": True,
        "remaining_label_bank_complete_at_freeze": False,
        "candidate_outcomes_used_for_implementation": False,
        "source_hashes": {name: sha256(ROOT / name) for name in FILES},
        "supersedes": "frozen_implementation_manifest_v6.json",
        "post_v6_change_scope": "execution moves outside the filesystem sandbox because sandboxed MATLAB fails ARM CPU-feature detection; the four shards remain sequential and all scientific definitions are unchanged",
        "note": "Any post-freeze source change must be versioned and reported before the corresponding result is interpreted as prespecified.",
    }
    OUTFILE.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
