#!/usr/bin/env python3
"""Paired halo comparison of fixed b=0.10 and 30 pkpc aperture."""
from pathlib import Path
import csv
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
SUMMARY = ROOT / "supplement"
INPUT = SUMMARY / "camels_host_metrics.csv"
OUTPUT = SUMMARY / "camels_b010_vs_aperture.csv"
SEED = 20261030
NBOOT = 5000
METRICS = [
    "median_iou", "median_completeness", "median_contamination_fraction_outside_contour",
    "median_log10_luminosity_ratio", "median_abs_log10_luminosity_error",
    "median_log10_R50_ratio", "median_abs_log10_R50_error",
]

frame = pd.read_csv(INPUT)
fof = frame[(frame.method == "FoF") & (frame.b_star.astype(str) == "0.10")]
ap = frame[frame.method == "aperture30pkpc"]
keys = ["cv_box", "fof_group_index", "band"]
paired = fof.merge(ap, on=keys, suffixes=("_b010", "_aperture"), validate="one_to_one")
if len(paired) != 98:
    raise RuntimeError(f"Expected 49 paired hosts in each band, found {len(paired)} rows")
if not (paired.n_failed_projections_b010.eq(0).all() and paired.n_failed_projections_aperture.eq(0).all()):
    raise RuntimeError("Unexpected failed projection in paired comparison")

rng = np.random.default_rng(SEED)
rows = []
box_ids = sorted(paired.cv_box.unique())
for band in ("V", "g"):
    subset = paired[paired.band == band].copy().sort_values(keys[:2])
    if len(subset) != 49:
        raise RuntimeError(f"Incomplete halo set in band {band}")
    for metric in METRICS:
        delta = (subset[f"{metric}_b010"] - subset[f"{metric}_aperture"]).to_numpy(float)
        host_draws = rng.choice(delta, size=(NBOOT, len(delta)), replace=True).mean(axis=1)
        box_draws = np.empty(NBOOT)
        by_box = {b: subset[subset.cv_box == b][metric + "_b010"].to_numpy(float)
                  - subset[subset.cv_box == b][metric + "_aperture"].to_numpy(float)
                  for b in box_ids}
        for i in range(NBOOT):
            selected = rng.choice(box_ids, size=len(box_ids), replace=True)
            draw = np.concatenate([rng.choice(by_box[b], size=len(by_box[b]), replace=True)
                                   for b in selected])
            box_draws[i] = draw.mean()
        rows.append({
            "band": band, "metric": metric, "n_paired_hosts": int(len(delta)),
            "mean_paired_difference_b010_minus_aperture": float(delta.mean()),
            "host_bootstrap_95_low": float(np.quantile(host_draws, .025)),
            "host_bootstrap_95_high": float(np.quantile(host_draws, .975)),
            "box_then_host_bootstrap_95_low": float(np.quantile(box_draws, .025)),
            "box_then_host_bootstrap_95_high": float(np.quantile(box_draws, .975)),
        })

pd.DataFrame(rows).to_csv(OUTPUT, index=False, quoting=csv.QUOTE_MINIMAL, float_format="%.8g")
print(pd.DataFrame(rows).to_string(index=False))
