#!/usr/bin/env python3
"""Consolidate pair yield, compute cost, and exact-edge recovery metrics."""

import json
import os
from datetime import datetime, timezone
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/private/tmp/brusselator_matplotlib")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parent
PAIR_EVALUATION = ROOT / "experiment_outputs" / "ml_BC_pair_acquisition_evaluation_v1"
PAIR_LABELING = ROOT / "experiment_outputs" / "ml_BC_pair_acquisition_labeling_v1"
PAIR_RECOVERY = ROOT / "experiment_outputs" / "ml_BC_pair_acquisition_edge_recovery_v1"
REFERENCE_RECOVERY = ROOT / "experiment_outputs" / "ml_BC_edge_recovery_v1"
OUTDIR = ROOT / "experiment_outputs" / "ml_BC_pair_acquisition_benchmark_v1"


def metadata_runtime(path):
    data = json.loads(path.read_text(encoding="utf-8"))
    return float(data["metadata"]["runtime_seconds"])


def main():
    if OUTDIR.exists():
        raise FileExistsError(f"Refusing to overwrite benchmark: {OUTDIR}")
    pair_outcomes = pd.read_csv(PAIR_EVALUATION / "pair_endpoint_outcomes.csv")
    pair_summary = pd.read_csv(PAIR_RECOVERY / "edge_recovery_summary.csv")
    reference = pd.read_csv(REFERENCE_RECOVERY / "edge_recovery_summary.csv")
    endpoint_labels = pd.read_csv(PAIR_LABELING / "pair_acquisition_endpoint_labels.csv")
    bisections = pd.read_csv(PAIR_RECOVERY / "bisection_cases.csv")

    methods = pd.DataFrame(
        {
            "method": ["pointwise_active", "space_filling", "pair_aware"],
            "frozen_pair_budget": [4, 4, 4],
            "new_endpoint_PDE_runs": [8, 8, 8],
            "verified_BC_brackets": [
                int(reference.loc[reference.model == "active", "verified_BC_bracket_found"].iloc[0]),
                int(reference.loc[reference.model == "baseline", "verified_BC_bracket_found"].iloc[0]),
                int(pair_outcomes.verified_BC_bracket.sum()),
            ],
            "exact_EBC_recoveries": [
                int(reference.loc[reference.model == "active", "recovered_exact_edge_neighborhood"].iloc[0]),
                int(reference.loc[reference.model == "baseline", "recovered_exact_edge_neighborhood"].iloc[0]),
                int(pair_summary.recovered_exact_edge_neighborhood.sum()),
            ],
        }
    )
    methods["bracket_yield_per_endpoint_run"] = methods.verified_BC_brackets / methods.new_endpoint_PDE_runs
    methods["end_to_end_recovery_per_endpoint_run"] = methods.exact_EBC_recoveries / methods.new_endpoint_PDE_runs
    methods["recovery_fraction_given_bracket"] = np.where(
        methods.verified_BC_brackets > 0,
        methods.exact_EBC_recoveries / methods.verified_BC_brackets,
        np.nan,
    )

    successful = pair_summary[pair_summary.recovered_exact_edge_neighborhood.astype(bool)].copy()
    endpoint_runtime = float(endpoint_labels.runtime_seconds.sum())
    bisection_runtime = float(bisections.runtime_seconds.sum())
    downstream_runtime = metadata_runtime(PAIR_RECOVERY / "configuration_and_metadata.json")
    total_runtime = endpoint_runtime + downstream_runtime
    summary = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "status": "completed retrospective B-C method-development benchmark",
        "pair_acquisition": {
            "frozen_pairs": 4,
            "new_endpoint_PDE_runs": 8,
            "verified_BC_brackets": int(pair_outcomes.verified_BC_bracket.sum()),
            "bracket_yield_per_endpoint_run": float(pair_outcomes.verified_BC_bracket.sum() / 8),
        },
        "dynamical_usefulness": {
            "valid_brackets_tested": int(pair_summary.verified_BC_bracket_found.sum()),
            "exact_EBC_recoveries": int(pair_summary.recovered_exact_edge_neighborhood.sum()),
            "recovery_fraction_given_valid_bracket": float(
                pair_summary.recovered_exact_edge_neighborhood.sum()
                / pair_summary.verified_BC_bracket_found.sum()
            ),
            "minimum_edge_orbit_distance_range": [
                float(successful.minimum_edge_orbit_distance.min()),
                float(successful.minimum_edge_orbit_distance.max()),
            ],
            "edge_shadow_duration_range": [
                float(successful.edge_shadow_duration.min()),
                float(successful.edge_shadow_duration.max()),
            ],
            "direct_period_relative_error_range": [
                float(successful.period_relative_error.min()),
                float(successful.period_relative_error.max()),
            ],
            "exact_edge_period": float(successful.exact_edge_period.iloc[0]),
        },
        "observed_compute_cost_seconds": {
            "eight_endpoint_labels": endpoint_runtime,
            "thirty_bisection_classifications": bisection_runtime,
            "complete_downstream_recovery_batch": downstream_runtime,
            "endpoint_plus_complete_downstream": total_runtime,
        },
        "experimental_separation": {
            "pair_selection_used_endpoint_outcomes": False,
            "pair_selection_used_exact_edge_orbit": False,
            "pair_selection_used_held_out_CA_region": False,
            "fixed_test_outcomes_used_for_selection": False,
            "fixed_test_state_indices_used_for_exclusion_only": True,
            "interpretation": (
                "strong method-development evidence, not a held-out confirmation; "
                "four proposals are too few for inferential claims"
            ),
        },
    }

    OUTDIR.mkdir(parents=True)
    methods.to_csv(OUTDIR / "method_comparison.csv", index=False)
    successful.to_csv(OUTDIR / "pair_aware_exact_edge_recoveries.csv", index=False)
    with (OUTDIR / "benchmark_summary.json").open("w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2)

    figure, axes = plt.subplots(1, 3, figsize=(13, 4.2), constrained_layout=True)
    x = np.arange(len(methods))
    width = 0.36
    axes[0].bar(x - width / 2, methods.verified_BC_brackets, width, label="verified bracket")
    axes[0].bar(x + width / 2, methods.exact_EBC_recoveries, width, label="exact-edge recovery")
    axes[0].set_xticks(x, ["point-active", "space-fill", "pair-aware"])
    axes[0].set_ylim(0, 4.5)
    axes[0].set_ylabel("count out of four frozen pairs")
    axes[0].legend()
    axes[0].set_title("Downstream usefulness")

    axes[1].semilogy(successful.model, successful.minimum_edge_orbit_distance, "o", markersize=8)
    axes[1].axhline(2e-2, color="black", linestyle="--", label="recovery gate")
    axes[1].set_ylabel("minimum distance to exact $E_{BC}$")
    axes[1].set_title("Exact-orbit proximity")
    axes[1].legend()

    axes[2].bar(successful.model, successful.edge_shadow_duration)
    axes[2].set_ylabel("time in exact-edge neighborhood")
    axes[2].set_title("Edge shadowing")
    figure.suptitle("Joint-posterior pair acquisition: B--C method-development benchmark")
    figure.savefig(OUTDIR / "pair_acquisition_dynamical_benchmark.png", dpi=300)
    plt.close(figure)

    print(json.dumps(summary, indent=2))
    print(methods.to_string(index=False))


if __name__ == "__main__":
    main()
