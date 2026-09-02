#!/usr/bin/env python3
"""Deduplicate conditionally selected physical pairs before edge tracking.

This is a compute-only transformation: every method/replicate keeps its frozen
highest-ranked-valid selection, but an identical physical pair is integrated
once and its dynamical result is attributed back to every selecting group.
"""

from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parent
EVALUATION_DIR = ROOT / "experiment_outputs" / "ml_BC_matched_acquisition_evaluation_v1"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    source = EVALUATION_DIR / "selected_edge_recovery_pairs.csv"
    unique_path = EVALUATION_DIR / "unique_edge_recovery_pairs.csv"
    attribution_path = EVALUATION_DIR / "edge_recovery_attribution.csv"
    manifest_path = EVALUATION_DIR / "unique_edge_recovery_manifest.json"
    if unique_path.exists() or attribution_path.exists() or manifest_path.exists():
        raise FileExistsError("Refusing to overwrite unique edge-recovery manifests")
    selected = pd.read_csv(source)
    selected["canonical_first_state_index"] = np.minimum(
        selected.first_state_index, selected.second_state_index
    ).astype(int)
    selected["canonical_second_state_index"] = np.maximum(
        selected.first_state_index, selected.second_state_index
    ).astype(int)
    selected["edge_run_id"] = [
        f"edge_pair_{left:03d}_{right:03d}"
        for left, right in zip(
            selected.canonical_first_state_index, selected.canonical_second_state_index
        )
    ]
    unique = selected.drop_duplicates(
        ["canonical_first_state_index", "canonical_second_state_index"], keep="first"
    ).copy()
    unique.to_csv(unique_path, index=False)
    selected.to_csv(attribution_path, index=False)
    payload = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "selection_count": len(selected),
        "unique_physical_pair_count": len(unique),
        "duplicate_integration_count_avoided": len(selected) - len(unique),
        "scientific_selection_changed": False,
        "selection_rule": "retain every frozen method/replicate selection; integrate each identical physical endpoint pair once",
        "selected_edge_recovery_pairs_sha256": sha256(source),
        "unique_edge_recovery_pairs_sha256": sha256(unique_path),
        "edge_recovery_attribution_sha256": sha256(attribution_path),
    }
    manifest_path.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
