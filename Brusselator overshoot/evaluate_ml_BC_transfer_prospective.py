#!/usr/bin/env python3
"""Open frozen selected-endpoint outcomes and evaluate transfer success."""

from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parent
CONFIG = ROOT / "ml_bc_transfer_prospective_config.json"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    config = json.loads(CONFIG.read_text(encoding="utf-8"))
    outdir = ROOT / config["evaluation_directory"]
    if outdir.exists():
        raise FileExistsError(f"Refusing to overwrite transfer evaluation: {outdir}")
    proposal_dir = ROOT / config["proposal_directory"]
    manifest_path = proposal_dir / "frozen_transfer_proposal_manifest.json"
    proposal_path = proposal_dir / "frozen_transfer_pair_proposals.csv"
    endpoint_path = ROOT / config["endpoint_labels_csv"]
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if not manifest.get("proposal_batches_frozen") or manifest.get("candidate_outcomes_used_for_proposals"):
        raise ValueError("Transfer proposals were not frozen outcome-blind")
    if sha256(proposal_path) != manifest["proposal_table_sha256"]:
        raise ValueError("Transfer proposals changed before outcome opening")

    proposals = pd.read_csv(proposal_path)
    labels = pd.read_csv(endpoint_path, usecols=["state_index", "final_label", "raw_five_way_outcome", "runtime_seconds"])
    first = labels.rename(
        columns={
            "state_index": "first_state_index",
            "final_label": "first_outcome",
            "raw_five_way_outcome": "first_raw_outcome",
            "runtime_seconds": "first_runtime_seconds",
        }
    )
    second = labels.rename(
        columns={
            "state_index": "second_state_index",
            "final_label": "second_outcome",
            "raw_five_way_outcome": "second_raw_outcome",
            "runtime_seconds": "second_runtime_seconds",
        }
    )
    opened = proposals.merge(first, on="first_state_index", validate="many_to_one").merge(
        second, on="second_state_index", validate="many_to_one"
    )
    if len(opened) != len(proposals):
        raise ValueError("Not every frozen proposal endpoint was labeled")
    opened["verified_BC_bracket"] = [
        {left, right} == {"B", "C"} for left, right in zip(opened.first_outcome, opened.second_outcome)
    ]
    opened["contains_other_or_U"] = ~opened.first_outcome.isin(["B", "C"]) | ~opened.second_outcome.isin(["B", "C"])

    summary_rows = []
    selected_rows = []
    for method, frame in opened.groupby("method", sort=False):
        frame = frame.sort_values("proposal_rank")
        valid = frame[frame.verified_BC_bracket]
        summary_rows.append(
            {
                "method": method,
                "proposal_count": len(frame),
                "endpoint_budget": 2 * len(frame),
                "valid_BC_brackets": int(frame.verified_BC_bracket.sum()),
                "contains_other_or_U_pairs": int(frame.contains_other_or_U.sum()),
                "at_least_one_valid_bracket": not valid.empty,
                "first_valid_bracket_rank": int(valid.proposal_rank.min()) if not valid.empty else None,
                "endpoint_labels_to_first_bracket": 2 * int(valid.proposal_rank.min()) if not valid.empty else None,
            }
        )
        if not valid.empty:
            selected_rows.append(valid.iloc[0].to_dict())
    method_summary = pd.DataFrame(summary_rows)

    unique_rows: dict[tuple[int, int], dict] = {}
    for row in selected_rows:
        key = tuple(sorted([int(row["first_state_index"]), int(row["second_state_index"])]))
        if key not in unique_rows:
            item = dict(row)
            item["selecting_methods"] = str(row["method"])
            unique_rows[key] = item
        else:
            methods = set(unique_rows[key]["selecting_methods"].split(";")) | {str(row["method"])}
            unique_rows[key]["selecting_methods"] = ";".join(sorted(methods))
    edge_rows = []
    for position, item in enumerate(unique_rows.values(), start=1):
        edge_rows.append(
            {
                "edge_run_id": f"transfer_edge_pair_{position:02d}",
                "selecting_methods": item["selecting_methods"],
                "proposal_rank": int(item["proposal_rank"]),
                "first_state_index": int(item["first_state_index"]),
                "second_state_index": int(item["second_state_index"]),
                "first_probability_C": float(item["first_probability_C"]),
                "second_probability_C": float(item["second_probability_C"]),
                "first_outcome": item["first_outcome"],
                "second_outcome": item["second_outcome"],
                "normalized_separation": float(item["normalized_separation"]),
            }
        )
    edge_pairs = pd.DataFrame(edge_rows)

    outdir.mkdir(parents=True)
    opened_path = outdir / "opened_transfer_pair_proposals.csv"
    summary_path = outdir / "transfer_method_summary.csv"
    edge_path = outdir / "edge_recovery_pairs.csv"
    opened.to_csv(opened_path, index=False)
    method_summary.to_csv(summary_path, index=False)
    edge_pairs.to_csv(edge_path, index=False)
    summary = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "prospective_endpoint_outcomes_opened": True,
        "proposal_manifest_sha256": sha256(manifest_path),
        "endpoint_labels_sha256": sha256(endpoint_path),
        "opened_proposals_sha256": sha256(opened_path),
        "method_summary_sha256": sha256(summary_path),
        "edge_recovery_pairs_sha256": sha256(edge_path),
        "selected_endpoint_count": int(labels.state_index.nunique()),
        "selected_endpoint_runtime_seconds": float(labels.drop_duplicates("state_index").runtime_seconds.sum()),
        "third_or_unresolved_endpoint_count": int((~labels.final_label.isin(["B", "C"])).sum()),
        "methods_with_a_valid_bracket": int(method_summary.at_least_one_valid_bracket.sum()),
        "unique_pairs_selected_for_conditional_edge_recovery": len(edge_pairs),
        "exact_edge_orbit_used_for_acquisition": False,
        "interpretation": "single-family prospective transfer; descriptive rather than a repeated statistical benchmark",
    }
    (outdir / "transfer_evaluation_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(method_summary.to_string(index=False))
    print(json.dumps(summary, indent=2))
    print(f"Prospective transfer evaluation saved in: {outdir}")


if __name__ == "__main__":
    main()
