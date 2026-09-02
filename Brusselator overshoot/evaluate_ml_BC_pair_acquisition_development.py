#!/usr/bin/env python3
"""Evaluate the prospectively frozen pair batch after endpoint labels exist."""

import hashlib
import json
import os
from datetime import datetime, timezone
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/private/tmp/brusselator_matplotlib")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd


ROOT = Path(__file__).resolve().parent
FROZEN = ROOT / "experiment_outputs" / "ml_BC_pair_acquisition_frozen_v1"
LABELS = ROOT / "experiment_outputs" / "ml_BC_pair_acquisition_labeling_v1"
REFERENCE = ROOT / "experiment_outputs" / "ml_BC_edge_recovery_v1"
OUTDIR = ROOT / "experiment_outputs" / "ml_BC_pair_acquisition_evaluation_v1"


def sha256_file(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main():
    if OUTDIR.exists():
        raise FileExistsError(f"Refusing to overwrite existing evaluation: {OUTDIR}")
    manifest = json.loads((FROZEN / "frozen_pair_manifest.json").read_text())
    if not manifest["pair_batch_frozen"] or manifest["endpoint_outcomes_used"]:
        raise RuntimeError("Pair acquisition was not frozen before endpoint labeling")
    batch_file = FROZEN / "frozen_pair_batch.csv"
    if sha256_file(batch_file) != manifest["frozen_pair_batch_hash"]:
        raise RuntimeError("Frozen pair batch changed after labeling")

    pairs = pd.read_csv(batch_file)
    labels = pd.read_csv(LABELS / "pair_acquisition_endpoint_labels.csv")
    by_index = labels.set_index("state_index")
    pairs["first_actual_label"] = pairs.first_state_index.map(by_index.final_label)
    pairs["second_actual_label"] = pairs.second_state_index.map(by_index.final_label)
    pairs["verified_BC_bracket"] = pairs.apply(
        lambda row: {row.first_actual_label, row.second_actual_label} == {"B", "C"}, axis=1
    )
    pairs["endpoint_PDE_runs"] = 2
    pairs["brackets_per_endpoint_run"] = pairs.verified_BC_bracket.astype(float) / 2.0

    reference = pd.read_csv(REFERENCE / "edge_recovery_summary.csv")
    reference_yield = dict(zip(reference.model, reference.verified_BC_bracket_found.astype(int)))
    valid_count = int(pairs.verified_BC_bracket.sum())
    summary = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "phase": manifest["phase"],
        "pair_batch_was_frozen_before_endpoint_labels": True,
        "held_out_CA_region_used": False,
        "pair_count": int(len(pairs)),
        "new_endpoint_PDE_runs": int(2 * len(pairs)),
        "verified_BC_brackets": valid_count,
        "pair_success_fraction": valid_count / len(pairs),
        "verified_brackets_per_endpoint_PDE_run": valid_count / (2 * len(pairs)),
        "reference_four_pair_batches": {
            "pointwise_active_verified_brackets": int(reference_yield["active"]),
            "space_filling_verified_brackets": int(reference_yield["baseline"]),
            "pair_aware_verified_brackets": valid_count,
            "note": "Descriptive method-development comparison; n=4 pairs per method is not a significance test.",
        },
        "primary_next_metric": (
            "fraction of verified brackets whose unchanged bisection/candidate workflow "
            "recovers the independently Newton-converged E_BC orbit"
        ),
    }

    OUTDIR.mkdir(parents=True)
    pairs.to_csv(OUTDIR / "pair_endpoint_outcomes.csv", index=False)
    with (OUTDIR / "pair_acquisition_yield.json").open("w") as handle:
        json.dump(summary, handle, indent=2)

    names = ["Point-active", "Space-filling", "Pair-aware"]
    values = [reference_yield["active"], reference_yield["baseline"], valid_count]
    figure, axis = plt.subplots(figsize=(6.4, 4.2), constrained_layout=True)
    bars = axis.bar(names, values, color=["tab:green", "tab:cyan", "tab:orange"])
    axis.bar_label(bars, labels=[f"{value}/4" for value in values], padding=3)
    axis.set_ylim(0, 4.5)
    axis.set_ylabel("verified B--C endpoint pairs")
    axis.set_title("Frozen bracket yield (method-development plane)")
    figure.savefig(OUTDIR / "frozen_pair_yield.png", dpi=300)
    plt.close(figure)

    print(json.dumps(summary, indent=2))
    print(pairs[["pair_id", "first_actual_label", "second_actual_label", "verified_BC_bracket"]])


if __name__ == "__main__":
    main()
