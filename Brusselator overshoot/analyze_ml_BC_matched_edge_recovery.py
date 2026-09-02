#!/usr/bin/env python3
"""Aggregate conditional guarded-refinement and exact-edge recovery results."""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parent
EVALUATION_DIR = ROOT / "experiment_outputs" / "ml_BC_matched_acquisition_evaluation_v1"
RECOVERY_DIR = ROOT / "experiment_outputs" / "ml_BC_matched_edge_recovery_v1"
OUTDIR = ROOT / "experiment_outputs" / "ml_BC_matched_dynamical_benchmark_v1"


def main() -> None:
    if OUTDIR.exists():
        raise FileExistsError(f"Refusing to overwrite dynamical summary: {OUTDIR}")
    selected = pd.read_csv(EVALUATION_DIR / "edge_recovery_attribution.csv")
    recovery = pd.read_csv(RECOVERY_DIR / "edge_recovery_summary.csv")
    bisection = pd.read_csv(RECOVERY_DIR / "bisection_cases.csv")
    safe_rows = []
    for model, frame in bisection.groupby("model"):
        safe_rows.append(
            {
                "model": model,
                "guarded_bisection_step_count": len(frame),
                "multiclass_safe_bisection": len(frame) == 10 and set(frame.outcome) <= {"B", "C"},
            }
        )
    safe = pd.DataFrame(safe_rows)
    selected = selected.rename(columns={"edge_run_id": "model"})
    merged = selected.merge(recovery, on="model", how="left", validate="many_to_one").merge(
        safe, on="model", how="left", validate="many_to_one"
    )
    if merged.recovered_exact_edge_neighborhood.isna().any():
        raise ValueError("A conditionally selected edge test has no recovery result")
    merged["multiclass_safe_bisection"] = merged.multiclass_safe_bisection.fillna(False)

    rows = []
    for method, frame in merged.groupby("method"):
        rows.append(
            {
                "method": method,
                "conditionally_tested_valid_brackets": len(frame),
                "multiclass_safe_refinements": int(frame.multiclass_safe_bisection.sum()),
                "exact_EBC_recoveries": int(frame.recovered_exact_edge_neighborhood.sum()),
                "exact_EBC_recovery_fraction_given_tested_valid_bracket": float(
                    frame.recovered_exact_edge_neighborhood.mean()
                ),
                "median_minimum_edge_orbit_distance": float(frame.minimum_edge_orbit_distance.median()),
                "median_edge_shadow_duration": float(frame.edge_shadow_duration.median()),
            }
        )
    summary = pd.DataFrame(rows).sort_values("method")
    OUTDIR.mkdir(parents=True)
    merged.to_csv(OUTDIR / "conditional_edge_recovery_cases.csv", index=False)
    summary.to_csv(OUTDIR / "conditional_edge_recovery_summary.csv", index=False)
    payload = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "status": "conditional highest-ranked-valid-pair dynamical evaluation complete",
        "selection_rule": "highest-ranked independently verified B-C proposal per method and matched initial-design replicate",
        "method_results": summary.to_dict(orient="records"),
        "interpretation": "exact-edge recovery is conditional on bracket discovery and complements, but does not replace, the primary replicate-level bracket-yield comparison",
    }
    (OUTDIR / "matched_dynamical_benchmark_summary.json").write_text(
        json.dumps(payload, indent=2), encoding="utf-8"
    )
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
