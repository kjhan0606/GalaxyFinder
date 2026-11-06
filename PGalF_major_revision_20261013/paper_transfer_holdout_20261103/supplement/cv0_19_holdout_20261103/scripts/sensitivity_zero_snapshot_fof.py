#!/usr/bin/env python3
"""Post-hoc sensitivity excluding hosts with zero catalogue/raw-FOF ID overlap."""
from __future__ import annotations

import csv
from collections import defaultdict
from pathlib import Path

import numpy as np

from summarize_holdout import paired_contrasts, read_csv, sample_summary, write_csv


ROOT = Path(__file__).resolve().parent
PRIMARY = ROOT / "summary_cv019_20261103_v2_potential_centered"
OUT = ROOT / "summary_cv019_20261103_sensitivity_no_zero_fof_labels"


def main() -> None:
    if OUT.exists() and any(OUT.iterdir()):
        raise FileExistsError(f"Refusing to overwrite {OUT}")
    audit = read_csv(ROOT / "membership_audit.csv")
    excluded = {
        (int(row["cv_box"]), int(row["fof_group_index"]))
        for row in audit if float(row["group_id_jaccard"]) == 0.0
    }
    hosts = read_csv(PRIMARY / "host_metrics.csv")
    kept = [
        row for row in hosts
        if (int(row["cv_box"]), int(row["fof_group_index"])) not in excluded
    ]
    groups: dict[tuple[str, str], list[dict[str, object]]] = defaultdict(list)
    for row in kept:
        for key, value in row.items():
            if key in {"cv_box", "fof_group_index", "n_projections", "n_invalid_contours", "invalid_contour_any_view", "iou_below_0p90", "large_failure_host_median", "large_failure_any_orientation", "failure_invalid_contour", "failure_contamination", "failure_luminosity", "failure_R50"}:
                try:
                    row[key] = int(value)
                except (TypeError, ValueError):
                    pass
            elif key not in {"band", "method_key", "b_star"}:
                try:
                    row[key] = float(value)
                except (TypeError, ValueError):
                    pass
        groups[(str(row["band"]), str(row["method_key"]))].append(row)
    rng = np.random.default_rng(20261104)
    summaries = [sample_summary(rows, rng) for _, rows in sorted(groups.items())]
    contrasts = paired_contrasts(kept, seed=20261104)
    OUT.mkdir(parents=True, exist_ok=True)
    write_csv(OUT / "host_metrics.csv", kept)
    write_csv(OUT / "summary.csv", summaries)
    write_csv(OUT / "paired_contrasts.csv", contrasts)

    full = {
        (row["band"], row["method_key"]): row
        for row in read_csv(PRIMARY / "summary.csv")
    }
    sensitivity = {(row["band"], row["method_key"]): row for row in summaries}
    compare_fields = (
        "mean_host_median_iou",
        "mean_host_median_contamination",
        "p90_host_abs_log10_luminosity_error",
        "p90_host_abs_log10_R50_error",
        "fraction_hosts_large_failure_host_median",
        "fraction_hosts_iou_below_0p90",
    )
    comparison = []
    for key in sorted(sensitivity):
        row = {"band": key[0], "method_key": key[1], "n_hosts_full": 130, "n_hosts_sensitivity": len({(int(r["cv_box"]), int(r["fof_group_index"])) for r in kept if (r["band"], r["method_key"]) == key})}
        for field in compare_fields:
            a, b = float(full[key][field]), float(sensitivity[key][field])
            row[f"full_{field}"] = a
            row[f"sensitivity_{field}"] = b
            row[f"sensitivity_minus_full_{field}"] = b - a
        comparison.append(row)
    write_csv(OUT / "comparison_with_primary.csv", comparison)
    print(f"Excluded {len(excluded)} hosts with zero raw-FoF label overlap; retained {len({(int(r['cv_box']), int(r['fof_group_index'])) for r in kept})} hosts.")
    print(f"Wrote sensitivity summaries to {OUT}")


if __name__ == "__main__":
    main()
