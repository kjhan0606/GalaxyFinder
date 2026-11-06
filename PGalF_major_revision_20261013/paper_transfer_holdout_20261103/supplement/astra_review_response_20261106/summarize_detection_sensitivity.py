#!/usr/bin/env python3
"""Summarize Darwin overlap under a common central-pixel detection cut."""
from __future__ import annotations

import csv
from collections import defaultdict
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
PAPER = HERE.parents[1]
RESULTS = HERE / "results"
SCORES = RESULTS / "corrected_scores_per_view.tsv"
FLAGS = RESULTS / "detection_flags" / "central_detection_flags.tsv"
OUTPUT = RESULTS / "detection_sensitivity_summary.tsv"
BANDS = ("V", "g")
LOS = ("x", "y", "z")


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def quantiles(values):
    x = np.asarray(values, dtype=float)
    x = x[np.isfinite(x)]
    if not x.size:
        return (np.nan, np.nan, np.nan)
    return tuple(float(v) for v in np.quantile(x, (0.16, 0.50, 0.84)))


def main():
    scores = read_tsv(SCORES)
    flags = read_tsv(FLAGS)
    valid = defaultdict(dict)
    for row in flags:
        valid[(row["step"], int(row["nstar"]), row["band"])][row["los"]] = (
            int(row["valid_central_isophote"]) == 1
        )
    records = defaultdict(dict)
    for row in scores:
        key = (row["step"], int(row["nstar"]), row["band"], row["los"])
        records[key][round(float(row["b"]), 2)] = row

    output = []
    for band in BANDS:
        all_entries = sorted({(step, nstar) for step, nstar, b, los in records
                              if b == band})
        cohorts = {
            "all_75": all_entries,
            "central_pixel_detected_all_three_views": [
                key for key in all_entries
                if all(valid[(key[0], key[1], band)].get(los, False) for los in LOS)
            ],
        }
        for cohort, entries in cohorts.items():
            per_entry = defaultdict(list)
            for step, nstar in entries:
                for los in LOS:
                    view = records[(step, nstar, band, los)]
                    rows = list(view.values())
                    positional_best = min(rows, key=lambda r: (
                        -float(r["pos_iou"]) if np.isfinite(float(r["pos_iou"])) else np.inf,
                        float(r["b"])))
                    psf_best = min(rows, key=lambda r: (
                        -float(r["psf_iou"]) if np.isfinite(float(r["psf_iou"])) else np.inf,
                        float(r["b"])))
                    primary_b = round(float(view[0.10]["primary_coarse_best_b"]), 2)
                    at_primary = view[primary_b]
                    per_entry["b_primary"].append(float(view[0.10]["primary_coarse_best_b"]))
                    per_entry["b_pos_fine"].append(float(positional_best["b"]))
                    per_entry["b_psf_fine"].append(float(psf_best["b"]))
                    per_entry["pos_iou_at_primary"].append(float(at_primary["pos_iou"]))
                    per_entry["psf_iou_at_primary"].append(float(at_primary["psf_iou"]))
                    per_entry["pos_iou_at_010"].append(float(view[0.10]["pos_iou"]))
                    per_entry["pos_iou_at_020"].append(float(view[0.20]["pos_iou"]))
                    per_entry["psf_iou_at_010"].append(float(view[0.10]["psf_iou"]))
                    per_entry["psf_iou_at_020"].append(float(view[0.20]["psf_iou"]))
                    per_entry["psf_iou_at_psf_fine_best"].append(float(psf_best["psf_iou"]))
                # First reduce the three projections, then summarize halo entries.
            row = {"band": band, "cohort": cohort, "n_halo_outputs": len(entries)}
            for statistic, values in per_entry.items():
                host_values = [float(np.nanmedian(values[i:i+3]))
                               for i in range(0, len(values), 3)]
                q16, median, q84 = quantiles(host_values)
                row[f"{statistic}_p16"] = q16
                row[f"{statistic}_median"] = median
                row[f"{statistic}_p84"] = q84
            output.append(row)

    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    tmp = OUTPUT.with_suffix(OUTPUT.suffix + ".tmp")
    with tmp.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(output[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(output)
    tmp.replace(OUTPUT)
    for row in output:
        print(row["band"], row["cohort"], row["n_halo_outputs"],
              "b primary", row["b_primary_median"],
              "b positional fine", row["b_pos_fine_median"],
              "b PSF fine", row["b_psf_fine_median"],
              "PSF IoU at primary", row["psf_iou_at_primary_median"],
              "PSF IoU at PSF optimum", row["psf_iou_at_psf_fine_best_median"])
    print(OUTPUT)


if __name__ == "__main__":
    main()
