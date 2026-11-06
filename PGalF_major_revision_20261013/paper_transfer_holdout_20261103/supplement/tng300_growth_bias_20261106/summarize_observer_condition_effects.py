#!/usr/bin/env python3
"""Pair actual-redshift and fixed-redshift growth biases by descendant tree."""
from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
INPUT = HERE / "tree_growth_bias_per_tree.csv"
OUTPUT = HERE / "observer_condition_growth_effects.csv"
SEED = 20261106
N_BOOT = 20000


def main():
    if OUTPUT.exists():
        raise FileExistsError(f"Refusing to overwrite {OUTPUT.name}")
    rows = list(csv.DictReader(INPUT.open(newline="")))
    by_key = {}
    for row in rows:
        key = (int(row["sample_rank"]), row["band"], row["method"], row["condition"])
        by_key[key] = row
    ranks = sorted({int(row["sample_rank"]) for row in rows})
    rng = np.random.default_rng(SEED)
    output = []
    for band in ("V", "g"):
        methods = sorted({row["method"] for row in rows if row["band"] == band})
        for method in methods:
            for metric, field in (("stellar_mass_growth_bias", "mass_growth_bias_dex"),
                                  ("R50_growth_bias", "R50_growth_bias_dex")):
                paired = np.array([
                    float(by_key[(rank, band, method, "actual_observer_redshift")][field])
                    - float(by_key[(rank, band, method, "fixed_observer_z0p197")][field])
                    for rank in ranks
                ])
                samples = paired[rng.integers(0, len(paired), size=(N_BOOT, len(paired)))]
                mean_ci = np.quantile(np.mean(samples, axis=1), [0.025, 0.975])
                median_ci = np.quantile(np.median(samples, axis=1), [0.025, 0.975])
                output.append({
                    "band": band,
                    "method": method,
                    "metric": metric,
                    "contrast": "actual observer redshift minus fixed observer redshift 0.1973",
                    "n_independent_descendant_trees": len(ranks),
                    "mean_paired_change_dex": float(np.mean(paired)),
                    "mean_paired_change_ci95": json.dumps([float(v) for v in mean_ci]),
                    "median_paired_change_dex": float(np.median(paired)),
                    "median_paired_change_ci95": json.dumps([float(v) for v in median_ci]),
                })
    with OUTPUT.open("x", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(output[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(output)
    print(json.dumps({"output": OUTPUT.name, "n_trees": len(ranks),
                      "rows": len(output)}, indent=2))


if __name__ == "__main__":
    main()
