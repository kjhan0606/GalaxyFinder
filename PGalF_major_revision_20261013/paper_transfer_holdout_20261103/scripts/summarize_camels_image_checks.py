#!/usr/bin/env python3
"""Reproduce fixed-parameter image checks from the supplied projection table.

Run with numpy and pandas, without snapshot files, yt, or FSPS. Hosts, rather
than views or noise realisations, are the units of statistical inference.
"""
from __future__ import annotations
import argparse
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
NBOOT = 5000
HOST_KEYS = ["cv_box", "fof_group_index"]
METRICS = ("iou", "completeness", "contamination", "logL", "abs_logL",
           "logR50", "abs_logR50", "failure_fraction")
CONDITIONS = {"baseline", "limit_brighter", "limit_fainter", "pixel_coarse",
              "no_psf_wing", "sky_residual_29", "noise_5sigma", "noise_3sigma",
              "noise_1sigma"}


def resampling_weights(hosts):
    rng = np.random.default_rng(20261030)
    n = len(hosts)
    plain = rng.multinomial(n, np.full(n, 1/n), size=NBOOT).astype(float)/n
    boxes = sorted(hosts.cv_box.unique())
    by_box = {b: np.flatnonzero(hosts.cv_box.to_numpy() == b) for b in boxes}
    nested = np.zeros((NBOOT, n))
    for k in range(NBOOT):
        for b in rng.choice(boxes, len(boxes), replace=True):
            indices = by_box[b]
            np.add.at(nested[k], rng.choice(indices, len(indices), replace=True), 1)
    nested /= nested.sum(axis=1, keepdims=True)
    return plain, nested


def paired(left, right, weights, metadata):
    merged = left.merge(right, on=HOST_KEYS, validate="one_to_one",
                        suffixes=("_left", "_right")).sort_values(HOST_KEYS)
    if len(merged) != 49:
        raise ValueError("Paired contrast lost hosts")
    rows = []
    for metric in METRICS:
        values = (merged[metric+"_left"]-merged[metric+"_right"]).to_numpy(float)
        good = np.isfinite(values)
        row = dict(metadata, metric=metric, n_hosts=49,
                   n_hosts_with_valid_metric=int(good.sum()),
                   mean_paired_difference=float(np.mean(values[good])) if good.any() else np.nan)
        for label, w in zip(("host", "box_then_host"), weights):
            denominator = w @ good.astype(float)
            means = np.divide(w @ np.nan_to_num(values, nan=0), denominator,
                              out=np.full(NBOOT, np.nan), where=denominator > 0)
            finite = np.isfinite(means)
            row[label+"_bootstrap95_low"] = float(np.quantile(means[finite], .025)) if finite.any() else np.nan
            row[label+"_bootstrap95_high"] = float(np.quantile(means[finite], .975)) if finite.any() else np.nan
        rows.append(row)
    return rows


def summarize(raw):
    raw = raw.copy()
    raw.b_star = raw.b_star.astype(str)
    keys = HOST_KEYS+["condition", "noise_realisation", "line_of_sight", "band", "method", "b_star"]
    if len(raw) != 13230 or set(raw.condition) != CONDITIONS or raw.duplicated(keys).any():
        raise ValueError("Incomplete or duplicated nine-condition input")
    failed = raw.target_failed.astype(bool)
    raw["failure_fraction"] = failed.astype(float)
    raw["contamination"] = (1-raw.precision).where(~failed)
    # Raw metrics explicitly encode zero IoU/completeness for a failed target.
    if not (raw.loc[failed, "iou"].eq(0).all() and raw.loc[failed, "completeness"].eq(0).all()):
        raise ValueError("Failed-target convention changed")
    for source, target in [("L", "logL"), ("R50", "logR50")]:
        num, den = ((raw.L_selected, raw.L_contour) if source == "L" else
                    (raw.R50_selected_kpc, raw.R50_contour_kpc))
        valid = (~failed) & num.gt(0) & den.gt(0) & np.isfinite(num) & np.isfinite(den)
        raw[target] = np.log10((num/den).where(valid))
        raw["abs_"+target] = raw[target].abs()
    raw["image_area_arcsec2"] = raw.target_area_arcsec2.where(~failed)
    real_keys = HOST_KEYS+["condition", "noise_realisation", "band", "method", "b_star"]
    if not raw.groupby(real_keys).size().eq(3).all():
        raise ValueError("A host-realisation lacks three projections")
    aggregations = {k: "median" for k in METRICS if k != "failure_fraction"}
    aggregations.update(failure_fraction="mean", image_area_arcsec2="median")
    real = raw.groupby(real_keys, as_index=False).agg(aggregations)
    host_keys = HOST_KEYS+["condition", "band", "method", "b_star"]
    host = real.groupby(host_keys, as_index=False)[list(METRICS)+["image_area_arcsec2"]].mean()
    host["n_noise_realisations"] = real.groupby(host_keys).size().to_numpy()
    expected = np.where(host.condition.str.startswith("noise_"), 3, 1)
    if not np.array_equal(host.n_noise_realisations.to_numpy(), expected):
        raise ValueError("Noise repetition counts differ from the design")
    summary = []
    for (condition, band, method, b), g in host.groupby(["condition", "band", "method", "b_star"]):
        if len(g) != 49:
            raise ValueError("Condition lacks 49 host entries")
        summary.append(dict(
            condition=condition, band=band, method=method, b_star=b, n_hosts=49,
            mean_host_iou=g.iou.mean(), median_host_iou=g.iou.median(),
            fraction_hosts_iou_below_0p90=g.iou.lt(.9).mean(),
            fraction_hosts_any_failed_view=g.failure_fraction.gt(0).mean(),
            mean_failed_projection_fraction=g.failure_fraction.mean(),
            median_host_completeness=g.completeness.median(),
            median_host_contamination=g.contamination.median(),
            mean_signed_logL_bias=g.logL.mean(), p90_absolute_logL_error=g.abs_logL.quantile(.9),
            mean_signed_logR50_bias=g.logR50.mean(), p90_absolute_logR50_error=g.abs_logR50.quantile(.9),
            n_hosts_valid_luminosity_error=int(g.abs_logL.notna().sum()),
            n_hosts_valid_radius_error=int(g.abs_logR50.notna().sum())))
    hosts = host[HOST_KEYS].drop_duplicates().sort_values(HOST_KEYS).reset_index(drop=True)
    if len(hosts) != 49 or hosts.cv_box.nunique() != 7:
        raise ValueError("Host manifest mismatch")
    weights = resampling_weights(hosts)
    changes, methods = [], []
    for (condition, band, method, b), g in host.groupby(["condition", "band", "method", "b_star"]):
        if condition == "baseline":
            continue
        baseline = host[host.condition.eq("baseline") & host.band.eq(band) & host.method.eq(method) & host.b_star.eq(b)]
        changes.extend(paired(g, baseline, weights, dict(condition=condition, band=band, method=method,
                            b_star=b, difference="condition minus baseline")))
    for (condition, band), g in host.groupby(["condition", "band"]):
        fixed = g[g.method.eq("FoF") & g.b_star.eq("0.1")]
        for method, b in [("FoF", "0.2"), ("aperture30pkpc", "aperture")]:
            other = g[g.method.eq(method) & g.b_star.eq(b)]
            methods.extend(paired(fixed, other, weights, dict(condition=condition, band=band,
                                  baseline_method=method, baseline_b_star=b,
                                  difference="fixed b=0.10 minus comparison method")))
    return {"host_realisation_metrics.csv": real, "host_metrics.csv": host,
            "sensitivity_summary.csv": pd.DataFrame(summary),
            "paired_condition_changes.csv": pd.DataFrame(changes),
            "paired_method_contrasts.csv": pd.DataFrame(methods)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=ROOT/"supplement/camels_image_checks/projection_metrics.csv.gz")
    parser.add_argument("--output", type=Path, required=True, help="New output directory")
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError("Use a new directory to preserve archived outputs")
    products = summarize(pd.read_csv(args.input))
    args.output.mkdir(parents=True)
    for name, frame in products.items():
        frame.to_csv(args.output/name, index=False)
    print(f"Reproduced {len(products)} tables from 49 hosts and nine conditions")


if __name__ == "__main__":
    main()
