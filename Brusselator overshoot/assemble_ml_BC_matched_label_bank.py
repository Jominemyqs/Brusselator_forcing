#!/usr/bin/env python3
"""Assemble the complete frozen-plane label bank with provenance checks."""

from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parent
CONFIG_FILE = ROOT / "ml_bc_matched_benchmark_protocol_config.json"
PROTOCOL_DIR = ROOT / "experiment_outputs" / "ml_BC_matched_benchmark_protocol_v1"
SHARD_PLAN_DIR = ROOT / "experiment_outputs" / "ml_BC_matched_label_bank_shards_v1"
MISSING_LABEL_SOURCES = [
    ROOT / "experiment_outputs" / "ml_BC_matched_label_bank_missing_v1" / "labeling_progress.csv",
] + [
    ROOT
    / "experiment_outputs"
    / f"ml_BC_matched_label_bank_shard{shard:02d}_v1"
    / f"matched_missing_labels_shard{shard:02d}.csv"
    for shard in range(1, 5)
]
OUTDIR = ROOT / "experiment_outputs" / "ml_BC_matched_label_bank_v1"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    if OUTDIR.exists():
        raise FileExistsError(f"Refusing to overwrite complete label bank: {OUTDIR}")
    config = json.loads(CONFIG_FILE.read_text(encoding="utf-8"))
    protocol_manifest_path = PROTOCOL_DIR / "frozen_protocol_manifest.json"
    protocol = json.loads(protocol_manifest_path.read_text(encoding="utf-8"))
    if not protocol.get("protocol_frozen") or protocol.get("candidate_outcomes_used"):
        raise ValueError("The matched protocol was not frozen outcome-blind")
    if sha256(PROTOCOL_DIR / "matched_initial_designs.csv") != protocol["matched_initial_designs_sha256"]:
        raise ValueError("Matched initial designs changed after protocol freeze")
    if sha256(PROTOCOL_DIR / "missing_label_manifest.csv") != protocol["missing_label_manifest_sha256"]:
        raise ValueError("Missing-label manifest changed after protocol freeze")
    shard_plan_path = SHARD_PLAN_DIR / "frozen_shard_plan.json"
    shard_plan = json.loads(shard_plan_path.read_text(encoding="utf-8"))
    if not shard_plan.get("shards_frozen") or shard_plan.get("outcomes_used_for_sharding"):
        raise ValueError("The label-bank shard plan was not frozen outcome-blind")
    if not all(path.exists() for path in MISSING_LABEL_SOURCES):
        raise FileNotFoundError("At least one checkpoint-prefix or shard label source is incomplete")

    sources = [ROOT / item for item in config["existing_label_sources"]] + MISSING_LABEL_SOURCES
    frames = []
    source_rows = []
    retained = [
        "state_index",
        "sample_id",
        "normalized_alpha_1",
        "normalized_alpha_2",
        "final_label",
        "label_status",
        "raw_five_way_outcome",
        "tested_duration",
        "decision_time",
        "runtime_seconds",
        "raw_file",
    ]
    for path in sources:
        frame = pd.read_csv(path)
        missing_columns = set(retained) - set(frame.columns)
        if missing_columns:
            raise ValueError(f"{path} is missing columns {sorted(missing_columns)}")
        item = frame[retained].copy()
        item["label_source"] = str(path.relative_to(ROOT))
        frames.append(item)
        source_rows.append(
            {
                "label_source": str(path.relative_to(ROOT)),
                "sha256": sha256(path),
                "row_count": len(item),
                "B_count": int((item.final_label == "B").sum()),
                "C_count": int((item.final_label == "C").sum()),
                "other_count": int((item.final_label == "other").sum()),
                "U_count": int((item.final_label == "U").sum()),
                "runtime_seconds": float(item.runtime_seconds.sum()),
            }
        )
    combined = pd.concat(frames, ignore_index=True)
    duplicate = combined[combined.duplicated("state_index", keep=False)]
    for state_index, group in duplicate.groupby("state_index"):
        if group.final_label.nunique() != 1:
            raise ValueError(f"Conflicting duplicate outcome for state {state_index}")
    combined = combined.sort_values("state_index").drop_duplicates("state_index", keep="first")

    geometry_path = ROOT / config["geometry_manifest_csv"]
    geometry = pd.read_csv(geometry_path)
    expected = int(config["geometry_gate"]["expected_state_count"])
    if len(combined) != expected or set(combined.state_index) != set(geometry.state_index):
        raise ValueError("Complete bank does not cover the frozen 17-by-17 plane")
    bank = geometry.merge(
        combined.drop(columns=["sample_id", "normalized_alpha_1", "normalized_alpha_2"]),
        on="state_index",
        how="left",
        validate="one_to_one",
    ).sort_values("state_index")

    OUTDIR.mkdir(parents=True)
    bank_path = OUTDIR / "complete_label_bank.csv"
    audit_path = OUTDIR / "label_source_audit.csv"
    bank.to_csv(bank_path, index=False)
    pd.DataFrame(source_rows).to_csv(audit_path, index=False)
    counts = bank.final_label.value_counts().to_dict()
    manifest = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "status": "complete",
        "state_count": len(bank),
        "B_count": int(counts.get("B", 0)),
        "C_count": int(counts.get("C", 0)),
        "other_count": int(counts.get("other", 0)),
        "U_count": int(counts.get("U", 0)),
        "binary_BC_plane_gate_passed": set(bank.final_label) <= {"B", "C"},
        "total_recorded_unique_PDE_runtime_seconds": float(bank.runtime_seconds.sum()),
        "complete_label_bank_sha256": sha256(bank_path),
        "label_source_audit_sha256": sha256(audit_path),
        "frozen_protocol_manifest_sha256": sha256(protocol_manifest_path),
        "frozen_shard_plan_sha256": sha256(shard_plan_path),
        "matched_initial_designs_sha256": protocol["matched_initial_designs_sha256"],
        "candidate_labels_are_for_masked_offline_evaluation_only": True,
        "interpretation": config["interpretation"],
    }
    (OUTDIR / "complete_label_bank_manifest.json").write_text(
        json.dumps(manifest, indent=2), encoding="utf-8"
    )
    print(json.dumps(manifest, indent=2))
    print(f"Complete matched label bank saved in: {OUTDIR}")


if __name__ == "__main__":
    main()
