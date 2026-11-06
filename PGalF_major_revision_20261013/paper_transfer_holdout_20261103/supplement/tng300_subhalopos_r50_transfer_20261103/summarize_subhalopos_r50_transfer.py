#!/usr/bin/env python3
"""Summarize TNG300 transfer metrics, including projected half-light radii."""
from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
SCORES = HERE / "host_metrics"
OUT = HERE / "host_level_summary"
METHODS = ("b020", "b010", "cal_primary", "cal_lowz_loo", "aperture30")
BANDS = ("V", "g")
CANDIDATES = ("b010", "cal_primary", "cal_lowz_loo")
COMPARATORS = ("b020", "b010", "aperture30")
LOG2 = float(np.log10(2.0))
N_BOOT = 10000
SEED = 20261103


def as_float(row, key):
    try:
        value = float(row[key])
    except (TypeError, ValueError, KeyError):
        return np.nan
    return value if np.isfinite(value) else np.nan


def read_view_rows():
    rows = []
    for rank in range(1, 41):
        path = SCORES / f"host_rank_{rank:02d}.csv"
        if not path.is_file():
            raise FileNotFoundError(f"Missing rank {rank} scores: {path}")
        with path.open(newline="") as f:
            rows.extend(csv.DictReader(f))
    expected = 40 * len(METHODS) * 3 * 2
    if len(rows) != expected:
        raise RuntimeError(f"Expected {expected} view-level rows, found {len(rows)}")
    return rows


def reduce_to_hosts(rows):
    grouped = {}
    for row in rows:
        key = (int(row["sample_rank"]), row["band"], row["method"])
        grouped.setdefault(key, []).append(row)
    expected_groups = 40 * len(METHODS) * len(BANDS)
    if len(grouped) != expected_groups:
        raise RuntimeError(f"Expected {expected_groups} host-method groups, found {len(grouped)}")
    hosts = []
    for (rank, band, method), views in sorted(grouped.items()):
        if len(views) != 3 or {r["los"] for r in views} != {"x", "y", "z"}:
            raise RuntimeError(f"Expected x/y/z views for rank={rank}, {band}, {method}")
        observed_edge = any(int(r["observed_edge"]) for r in views)
        observed_missing = any(int(r["observed_missing"]) for r in views)
        excluded = observed_edge or observed_missing

        def med(field):
            return float(np.nanmedian([as_float(r, field) for r in views]))

        def med_abs(field):
            return float(np.nanmedian(np.abs([as_float(r, field) for r in views])))

        log_l = med("log10_luminosity_ratio")
        log_r = med("log10_radius_ratio")
        log_r50 = med("log10_R50_ratio")
        contamination = med("contamination")
        lum_failure = not np.isfinite(log_l) or abs(log_l) > LOG2
        method_contour_failure = any(
            int(r["isophote_missing"]) or int(r["selected_edge"]) for r in views
        )
        radius_failure = (method_contour_failure or not np.isfinite(log_r)
                         or abs(log_r) > LOG2)
        r50_failure = (method_contour_failure or not np.isfinite(log_r50)
                       or abs(log_r50) > LOG2)
        contamination_failure = not np.isfinite(contamination) or contamination > 0.20
        hosts.append({
            "sample_rank": rank,
            "group_id": int(views[0]["group_id"]),
            "band": band,
            "method": method,
            "b": as_float(views[0], "b"),
            "observed_contour_excluded": int(excluded),
            "observed_edge_excluded": int(observed_edge),
            "observed_missing_excluded": int(observed_missing),
            "median_iou": med("iou"),
            "median_completeness": med("completeness"),
            "median_precision": med("precision"),
            "median_contamination": contamination,
            "median_log10_luminosity_ratio": log_l,
            "median_log10_radius_ratio": log_r,
            "median_log10_R50_ratio": log_r50,
            # Absolute errors must be computed per projection before reducing
            # the three repeated views to one host value. Taking abs(median(view
            # signed errors)) can cancel asymmetric view-level errors.
            "median_abs_log10_luminosity_ratio": med_abs("log10_luminosity_ratio"),
            "median_abs_log10_radius_ratio": med_abs("log10_radius_ratio"),
            "median_abs_log10_R50_ratio": med_abs("log10_R50_ratio"),
            "density_seed_offset_pkpc": med("density_seed_offset_pkpc"),
            "missing_or_edge_selected_views": sum(
                int(bool(int(r["isophote_missing"]) or int(r["selected_edge"]))) for r in views
            ),
            "luminosity_failure": int(lum_failure),
            "area_radius_failure": int(radius_failure),
            "R50_failure": int(r50_failure),
            "contamination_failure": int(contamination_failure),
            "large_failure": int(lum_failure or r50_failure or contamination_failure),
            "large_failure_area_radius": int(lum_failure or radius_failure or contamination_failure),
        })
    return hosts


def bootstrap_interval(values, statistic, rng):
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if values.size == 0:
        return np.nan, np.nan, np.nan
    point = float(statistic(values))
    draws = values[rng.integers(0, values.size, (N_BOOT, values.size))]
    boot = statistic(draws, axis=1)
    return point, float(np.quantile(boot, 0.025)), float(np.quantile(boot, 0.975))


def method_summaries(hosts):
    rng = np.random.default_rng(SEED)
    result = []
    for band in BANDS:
        for method in METHODS:
            all_hosts = [r for r in hosts if r["band"] == band and r["method"] == method]
            valid = [r for r in all_hosts if not r["observed_contour_excluded"]]
            if len(all_hosts) != 40:
                raise RuntimeError(f"Expected 40 hosts for {band}/{method}; found {len(all_hosts)}")
            row = {
                "band": band,
                "method": method,
                "b": all_hosts[0]["b"],
                "n_hosts_total": len(all_hosts),
                "n_hosts_valid": len(valid),
                "n_observed_edge_excluded": sum(r["observed_edge_excluded"] for r in all_hosts),
                "n_observed_missing_excluded": sum(r["observed_missing_excluded"] for r in all_hosts),
            }
            for field in (
                "median_iou", "median_completeness", "median_precision",
                "median_contamination", "median_log10_luminosity_ratio",
                "median_log10_radius_ratio", "median_log10_R50_ratio",
            ):
                values = np.asarray([r[field] for r in valid], dtype=float)
                values = values[np.isfinite(values)]
                row[f"mean_{field}"] = float(np.mean(values)) if values.size else np.nan
                row[f"median_{field}"] = float(np.median(values)) if values.size else np.nan
            for field in ("median_abs_log10_luminosity_ratio",
                          "median_abs_log10_radius_ratio", "median_abs_log10_R50_ratio"):
                values = np.asarray([r[field] for r in valid], dtype=float)
                values = values[np.isfinite(values)]
                row[f"p90_{field}"] = float(np.quantile(values, 0.90)) if values.size else np.nan
                point, lo, hi = bootstrap_interval(values, np.median, rng)
                row[f"median_{field}"] = point
                row[f"median_{field}_bootstrap_2p5"] = lo
                row[f"median_{field}_bootstrap_97p5"] = hi
            for field in ("median_iou", "median_contamination"):
                values = np.asarray([r[field] for r in valid], dtype=float)
                point, lo, hi = bootstrap_interval(values, np.median, rng)
                row[f"median_{field}_bootstrap_2p5"] = lo
                row[f"median_{field}_bootstrap_97p5"] = hi
            for field in ("large_failure", "large_failure_area_radius", "luminosity_failure",
                          "area_radius_failure", "R50_failure", "contamination_failure"):
                values = np.asarray([r[field] for r in valid], dtype=float)
                row[f"{field}_fraction"] = float(np.mean(values)) if values.size else np.nan
                row[f"{field}_count"] = int(np.sum(values)) if values.size else 0
            result.append(row)
    return result


def paired_contrasts(hosts):
    rng = np.random.default_rng(SEED + 1)
    out = []
    by_key = {(r["band"], r["method"], r["sample_rank"]): r for r in hosts}
    fields = (
        ("median_iou", "higher_is_better", "median"),
        ("median_completeness", "higher_is_better", "median"),
        ("median_precision", "higher_is_better", "median"),
        ("median_contamination", "lower_is_better", "median"),
        ("median_abs_log10_luminosity_ratio", "lower_is_better", "median"),
        ("median_abs_log10_radius_ratio", "lower_is_better", "median"),
        ("median_abs_log10_R50_ratio", "lower_is_better", "median"),
        ("large_failure", "lower_is_better", "mean"),
    )
    for band in BANDS:
        for comparator in COMPARATORS:
            for method in CANDIDATES:
                if method == comparator:
                    continue
                paired = [
                    (by_key[(band, comparator, rank)], by_key[(band, method, rank)])
                    for rank in range(1, 41)
                ]
                for field, direction, stat_name in fields:
                    diffs = np.asarray([
                        float(cand[field]) - float(ref[field])
                        for ref, cand in paired
                        if np.isfinite(float(ref[field])) and np.isfinite(float(cand[field]))
                    ])
                    if not diffs.size:
                        continue
                    stat = np.mean if stat_name == "mean" else np.median
                    draws = diffs[rng.integers(0, diffs.size, (N_BOOT, diffs.size))]
                    boot = stat(draws, axis=1)
                    out.append({
                        "band": band,
                        "method": method,
                        "comparator": comparator,
                        "metric": field,
                        "difference_candidate_minus_comparator": float(stat(diffs)),
                        "bootstrap_2p5": float(np.quantile(boot, 0.025)),
                        "bootstrap_97p5": float(np.quantile(boot, 0.975)),
                        "n_paired_hosts": int(diffs.size),
                        "interpretation": direction,
                        "bootstrap_unit": "host halo",
                        "bootstrap_resamples": N_BOOT,
                        "seed": SEED + 1,
                    })
    return out


def write_csv(path, rows):
    if not rows:
        raise RuntimeError(f"Refusing to write empty output: {path}")
    with path.open("x", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main():
    if OUT.exists():
        raise FileExistsError(f"Refusing to overwrite existing summary directory: {OUT}")
    OUT.mkdir(parents=True)
    view_rows = read_view_rows()
    hosts = reduce_to_hosts(view_rows)
    summaries = method_summaries(hosts)
    contrasts = paired_contrasts(hosts)
    write_csv(OUT / "host_level_metrics.csv", hosts)
    write_csv(OUT / "method_summary.csv", summaries)
    write_csv(OUT / "paired_method_contrasts.csv", contrasts)
    definitions = {
        "analysis_type": "SubhaloPos-centred follow-up with projected half-light-radius metrics after the original seed-centred TNG scores were inspected",
        "blindness_note": "TNG linking parameters were not tuned; this follow-up is not a fresh fully blinded preregistration because earlier seed-centred TNG scores were already open.",
        "calibration_parameters": {
            "darwin_primary_low_redshift": {"V": 0.12, "g": 0.10},
            "darwin_leave_one_host_out_low_redshift": {"V": 0.10, "g": 0.09},
            "fixed_controls": [0.10, 0.20],
        },
        "image_and_aperture_centre": "TNG SubhaloPos of the selected central subhalo (maximally bound particle position, ckpc/h); image centre is fixed across methods",
        "linking_root": "stellar-density peak, used only to choose the linked FoF component",
        "psf_fwhm_arcsec": 1.0,
        "statistical_unit": "FoF host halo",
        "view_reduction": "median over the three orthogonal projections within each host",
        "luminosity_error": "log10(L_selected / L_observed_inside_isophote)",
        "area_radius_error": "log10(R_selected_isophote / R_observed_isophote), with R=sqrt(area/pi)",
        "R50_error": "log10(R50_selected / R50_observed_contour); projected radii use full stellar flux at particle positions",
        "absolute_error_reduction": "median of absolute per-view log-ratios within each host; signed median ratios are reported separately",
        "contamination": "1 - light-weighted precision",
        "large_failure": "factor-of-two luminosity or R50 error, missing/edge selected contour, or contamination > 0.20",
        "observed_contour_exclusion": "host excluded for a band if any full-light isophote is absent or touches the image edge",
        "bootstrap": f"{N_BOOT} paired host resamples; seed {SEED + 1} for contrasts",
        "single_snapshot_limit": "one TNG300-1 snapshot and one independent simulation volume; no redshift-evolution or volume-to-volume claim",
    }
    with (OUT / "metric_definitions.json").open("x", encoding="utf-8") as f:
        json.dump(definitions, f, indent=2)
        f.write("\n")
    print(json.dumps({
        "view_rows": len(view_rows),
        "host_metric_rows": len(hosts),
        "method_summary_rows": len(summaries),
        "paired_contrast_rows": len(contrasts),
        "output": str(OUT),
    }, indent=2))


if __name__ == "__main__":
    main()
