#!/usr/bin/env python3
"""Freeze four disjoint execution shards for the remaining label bank."""

from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parent
PROTOCOL_DIR = ROOT / "experiment_outputs" / "ml_BC_matched_benchmark_protocol_v1"
PARTIAL_DIR = ROOT / "experiment_outputs" / "ml_BC_matched_label_bank_missing_v1"
OUTDIR = ROOT / "experiment_outputs" / "ml_BC_matched_label_bank_shards_v1"
SHARD_COUNT = 4


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    if OUTDIR.exists():
        raise FileExistsError(f"Refusing to overwrite frozen shard plan: {OUTDIR}")
    protocol = json.loads((PROTOCOL_DIR / "frozen_protocol_manifest.json").read_text())
    if not protocol.get("protocol_frozen") or protocol.get("candidate_outcomes_used"):
        raise ValueError("Parent matched protocol is not frozen outcome-blind")
    missing_path = PROTOCOL_DIR / "missing_label_manifest.csv"
    progress_path = PARTIAL_DIR / "labeling_progress.csv"
    missing = pd.read_csv(missing_path).sort_values("manifest_position")
    progress = pd.read_csv(progress_path)
    completed_indices = set(progress.state_index.astype(int))
    remaining = missing[~missing.state_index.isin(completed_indices)].copy().reset_index(drop=True)
    if len(missing) != protocol["missing_label_state_count"]:
        raise ValueError("Parent missing-label count changed")
    if len(completed_indices) + len(remaining) != len(missing):
        raise ValueError("Completed and remaining states do not partition the parent manifest")

    OUTDIR.mkdir(parents=True)
    completed = missing[missing.state_index.isin(completed_indices)].copy()
    completed.to_csv(OUTDIR / "completed_prefix_manifest.csv", index=False)
    shards = []
    for shard in range(1, SHARD_COUNT + 1):
        frame = remaining.iloc[(shard - 1) :: SHARD_COUNT].copy()
        frame.insert(0, "shard_position", range(1, len(frame) + 1))
        path = OUTDIR / f"shard_{shard:02d}_manifest.csv"
        frame.to_csv(path, index=False)
        shards.append(
            {
                "shard": shard,
                "path": str(path.relative_to(ROOT)),
                "state_count": len(frame),
                "sha256": sha256(path),
            }
        )
    combined = set()
    for item in shards:
        frame = pd.read_csv(ROOT / item["path"])
        values = set(frame.state_index.astype(int))
        if combined & values:
            raise ValueError("Frozen shards overlap")
        combined |= values
    if combined != set(remaining.state_index.astype(int)):
        raise ValueError("Frozen shards do not cover every remaining state")

    payload = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "shards_frozen": True,
        "shard_count": SHARD_COUNT,
        "parent_missing_label_state_count": len(missing),
        "completed_prefix_state_count": len(completed_indices),
        "remaining_state_count": len(remaining),
        "outcomes_used_for_sharding": False,
        "partition_rule": "round-robin by frozen parent manifest position after excluding checkpointed state indices",
        "parent_missing_label_manifest_sha256": sha256(missing_path),
        "completed_progress_sha256": sha256(progress_path),
        "shards": shards,
    }
    (OUTDIR / "frozen_shard_plan.json").write_text(json.dumps(payload, indent=2))
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
