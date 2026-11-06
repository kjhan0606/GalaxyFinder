#!/usr/bin/env python3
"""Summarize detection-qualified and full-cohort TNG growth measurements."""
from __future__ import annotations

import argparse
import csv
import json
from collections import defaultdict
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
DEFAULT_SCORES = HERE / "growth_scores_r1audit_20261106"
DEFAULT_OUTPUT = HERE / "growth_r1r4_audit_20261106"
SEED = 20261106
N_BOOT = 20000
SNAPS = (84, 67, 50)
ENDPOINTS = (50, 84)
LOS = ("x", "y", "z")
OBSERVERS = (0.19728418237601, 0.503047523244883, 0.99729422578194)
SNAP_TO_Z = {84: OBSERVERS[0], 67: OBSERVERS[1], 50: OBSERVERS[2]}
METHODS = (
    "isophotal_particle_selection",
    "b010_linked_group",
    "b012_linked_group",
    "b020_linked_group",
    "projected_aperture_30_pkpc",
    "Darwin_calibrated_sphere",
)
CONDITIONS = ("actual_observer_redshift", "fixed_observer_z0p197")
METRICS = {
    "stellar_mass": ("stellar_mass_msun", "central_bound_stellar_mass_msun"),
    "half_light_radius": ("projected_R50_pkpc", "central_bound_R50_pkpc"),
    "half_mass_radius": ("projected_R50_mass_pkpc", "central_bound_R50_mass_pkpc"),
}


def read_csv(path: Path):
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream))


def write_csv(path: Path, rows):
    if not rows:
        raise RuntimeError(f"No rows to write: {path}")
    with path.open("x", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def interval(values, rng):
    values = np.asarray(values, dtype=float)
    if values.size == 0:
        return None
    draws = values[rng.integers(0, values.size, size=(N_BOOT, values.size))]
    mean_draws = np.mean(draws, axis=1)
    median_draws = np.median(draws, axis=1)
    mae_draws = np.mean(np.abs(draws), axis=1)
    return {
        "n_trees_measurable": int(values.size),
        "mean_bias_dex": float(np.mean(values)),
        "mean_bias_ci95": np.quantile(mean_draws, [0.025, 0.975]).tolist(),
        "median_bias_dex": float(np.median(values)),
        "median_bias_ci95": np.quantile(median_draws, [0.025, 0.975]).tolist(),
        "mean_absolute_error_dex": float(np.mean(np.abs(values))),
        "mean_absolute_error_ci95": np.quantile(mae_draws, [0.025, 0.975]).tolist(),
        "p16_bias_dex": float(np.quantile(values, 0.16)),
        "p84_bias_dex": float(np.quantile(values, 0.84)),
    }


def as_float(row, key):
    value = row.get(key, "")
    try:
        return float(value)
    except (TypeError, ValueError):
        return float("nan")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--scores", type=Path, default=DEFAULT_SCORES)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)

    cohort = sorted(int(row["sample_rank"]) for row in read_csv(HERE / "lineage_cohort.csv")
                    if int(row["central_at_all_targets"]) == 1)
    indexed = {}
    detection_rows = []
    for rank in cohort:
        path = args.scores / f"tree_rank_{rank:02d}_growth_scores.csv"
        rows = read_csv(path)
        if len(rows) != 378:
            raise RuntimeError(f"Expected 378 rows in {path.name}; found {len(rows)}")
        for row in rows:
            snap = int(row["snapshot"])
            key = (rank, snap, row["band"], row["method"], row["line_of_sight"],
                   round(float(row["observer_redshift"]), 6))
            if key in indexed:
                raise RuntimeError(f"Duplicate score row: {key}")
            indexed[key] = row
            if row["method"] == "isophotal_particle_selection":
                detection_rows.append({
                    "sample_rank": rank,
                    "descendant_group_id": row["descendant_group_id"],
                    "snapshot": snap,
                    "physical_redshift": row["physical_redshift"],
                    "observer_redshift": row["observer_redshift"],
                    "band": row["band"],
                    "line_of_sight": row["line_of_sight"],
                    "status": row["isophote_status"],
                    "central_pixel_detected": row["central_pixel_detected"],
                    "nearest_bright_pixel_fallback_used": row["nearest_bright_pixel_fallback_used"],
                    "region_area_arcsec2": row["observed_isophote_area_arcsec2"],
                    "region_pixel_count": row["observed_isophote_pixel_count"],
                    "touches_image_edge": row["isophote_touches_image_edge"],
                    "centroid_offset_pixels": row["isophote_centroid_offset_pixels"],
                    "nearest_region_pixel_offset_pixels": row["nearest_isophote_pixel_offset_pixels"],
                    "selected_particle_count": row["selected_particle_count"],
                    "selected_stellar_mass_msun": row["stellar_mass_msun"],
                    "selected_half_light_radius_pkpc": row["projected_R50_pkpc"],
                })

    detection_summary = []
    groups = defaultdict(list)
    for row in detection_rows:
        groups[(int(row["snapshot"]), round(float(row["observer_redshift"]), 6),
                row["band"])].append(row)
    for (snap, observer_z, band), rows in sorted(groups.items()):
        n = len(rows)
        counts = {status: sum(row["status"] == status for row in rows)
                  for status in ("central_pixel_detected", "nearest_bright_pixel_fallback",
                                 "no_bright_region")}
        areas = [float(row["region_area_arcsec2"]) for row in rows
                 if float(row["region_area_arcsec2"]) > 0]
        detection_summary.append({
            "snapshot": snap,
            "physical_redshift": SNAP_TO_Z[snap],
            "observer_redshift": observer_z,
            "band": band,
            "n_tree_views": n,
            "n_central_pixel_detected": counts["central_pixel_detected"],
            "fraction_central_pixel_detected": counts["central_pixel_detected"] / n,
            "n_nearest_bright_pixel_fallback": counts["nearest_bright_pixel_fallback"],
            "fraction_fallback": counts["nearest_bright_pixel_fallback"] / n,
            "n_no_bright_region": counts["no_bright_region"],
            "fraction_no_bright_region": counts["no_bright_region"] / n,
            "n_region_touching_stamp_edge": sum(int(row["touches_image_edge"]) for row in rows),
            "median_detected_region_area_arcsec2": float(np.median(areas)) if areas else "",
        })
    write_csv(args.output / "r1_isophote_detection_by_epoch.csv", detection_summary)
    write_csv(args.output / "r1_isophote_detection_per_tree_view.csv", detection_rows)

    def get_row(rank, snap, band, method, los, condition):
        zobs = (SNAP_TO_Z[snap] if condition == "actual_observer_redshift" else OBSERVERS[0])
        key = (rank, snap, band, method, los, round(zobs, 6))
        try:
            return indexed[key]
        except KeyError as exc:
            raise RuntimeError(f"Missing score row {key}") from exc

    endpoint_valid_ranks = {}
    for band in ("V", "g"):
        for condition in CONDITIONS:
            valid = []
            for rank in cohort:
                detected = all(int(get_row(rank, snap, band, "isophotal_particle_selection",
                                           los, condition)["central_pixel_detected"]) == 1
                               for snap in ENDPOINTS for los in LOS)
                if detected:
                    valid.append(rank)
            endpoint_valid_ranks[(band, condition)] = valid

    rng = np.random.default_rng(SEED)
    growth_tree_rows = []
    growth_summary = []
    for band in ("V", "g"):
        for condition in CONDITIONS:
            detected_ranks = endpoint_valid_ranks[(band, condition)]
            for method in METHODS:
                view_biases = {metric: {} for metric in METRICS}
                status_by_rank = {}
                for rank in cohort:
                    n_views_with_any_failure = 0
                    per_metric = {metric: [] for metric in METRICS}
                    n_valid_views_by_metric = {metric: 0 for metric in METRICS}
                    for los in LOS:
                        start = get_row(rank, ENDPOINTS[0], band, method, los, condition)
                        end = get_row(rank, ENDPOINTS[1], band, method, los, condition)
                        detect_start = int(get_row(
                            rank, ENDPOINTS[0], band, "isophotal_particle_selection", los,
                            condition)["central_pixel_detected"])
                        detect_end = int(get_row(
                            rank, ENDPOINTS[1], band, "isophotal_particle_selection", los,
                            condition)["central_pixel_detected"])
                        if not (detect_start and detect_end):
                            n_views_with_any_failure += 1
                        for metric, (selected_key, reference_key) in METRICS.items():
                            vs, ve = as_float(start, selected_key), as_float(end, selected_key)
                            rs, re = as_float(start, reference_key), as_float(end, reference_key)
                            if min(vs, ve, rs, re) <= 0 or not np.all(np.isfinite([vs, ve, rs, re])):
                                continue
                            bias = float(np.log10(ve / vs) - np.log10(re / rs))
                            per_metric[metric].append(bias)
                            n_valid_views_by_metric[metric] += 1
                            view_biases[metric].setdefault(rank, []).append(bias)
                    status_by_rank[rank] = n_views_with_any_failure
                    row = {
                        "sample_rank": rank,
                        "band": band,
                        "condition": condition,
                        "method": method,
                        "n_endpoint_views_with_failed_central_detection": n_views_with_any_failure,
                        "detected_at_both_endpoints_in_all_three_views": int(rank in detected_ranks),
                    }
                    for metric in METRICS:
                        row[f"n_valid_view_pairs_{metric}"] = n_valid_views_by_metric[metric]
                        row[f"growth_bias_{metric}_dex"] = (
                            float(np.median(per_metric[metric])) if per_metric[metric] else "")
                    growth_tree_rows.append(row)

                for cohort_name, ranks in (("all_36_branches", cohort),
                                           ("central_pixel_detected_all_views_at_both_endpoints",
                                            detected_ranks)):
                    for metric in METRICS:
                        values = [float(np.median(view_biases[metric][rank]))
                                  for rank in ranks if rank in view_biases[metric]
                                  and view_biases[metric][rank]]
                        stats = interval(values, rng)
                        failures = sum(status_by_rank[rank] > 0 for rank in ranks)
                        summary_row = {
                            "band": band,
                            "observer_condition": condition,
                            "method": method,
                            "measurement": metric,
                            "cohort": cohort_name,
                            "n_trees_in_cohort": len(ranks),
                            "n_endpoint_detected_trees": len(detected_ranks),
                            "n_trees_with_any_endpoint_detection_failure": failures,
                            "bootstrap_unit": "descendant branch; median across available views per branch",
                        }
                        if stats:
                            summary_row.update(stats)
                        else:
                            summary_row.update({
                                "n_trees_measurable": 0, "mean_bias_dex": "",
                                "mean_bias_ci95": "", "median_bias_dex": "",
                                "median_bias_ci95": "", "mean_absolute_error_dex": "",
                                "mean_absolute_error_ci95": "", "p16_bias_dex": "",
                                "p84_bias_dex": "",
                            })
                        for key in ("mean_bias_ci95", "median_bias_ci95", "mean_absolute_error_ci95"):
                            if isinstance(summary_row.get(key), list):
                                summary_row[key] = json.dumps(summary_row[key])
                        growth_summary.append(summary_row)
    write_csv(args.output / "r1r3_growth_bias_per_tree.csv", growth_tree_rows)
    write_csv(args.output / "r1r3_growth_bias_summary.csv", growth_summary)

    mass_recovery_tree = []
    mass_recovery_summary = []
    for band in ("V", "g"):
        for condition in CONDITIONS:
            for snap in SNAPS:
                for method in METHODS:
                    tree_ratios = []
                    for rank in cohort:
                        view_ratios = []
                        for los in LOS:
                            row = get_row(rank, snap, band, method, los, condition)
                            selected_mass = as_float(row, "stellar_mass_msun")
                            central_mass = as_float(row, "central_bound_stellar_mass_msun")
                            ratio = selected_mass / central_mass if central_mass > 0 else float("nan")
                            if np.isfinite(ratio):
                                view_ratios.append(ratio)
                        tree_ratio = float(np.median(view_ratios)) if view_ratios else float("nan")
                        tree_ratios.append(tree_ratio)
                        mass_recovery_tree.append({
                            "sample_rank": rank, "band": band, "observer_condition": condition,
                            "snapshot": snap, "physical_redshift": SNAP_TO_Z[snap],
                            "method": method, "n_valid_projections": len(view_ratios),
                            "median_selected_to_central_mass_ratio": tree_ratio,
                        })
                    finite = np.asarray([v for v in tree_ratios if np.isfinite(v)], dtype=float)
                    mass_recovery_summary.append({
                        "band": band, "observer_condition": condition,
                        "snapshot": snap, "physical_redshift": SNAP_TO_Z[snap],
                        "method": method, "n_trees": int(finite.size),
                        "median_tree_level_mass_ratio": float(np.median(finite)) if finite.size else "",
                        "p16_tree_level_mass_ratio": float(np.quantile(finite, 0.16)) if finite.size else "",
                        "p84_tree_level_mass_ratio": float(np.quantile(finite, 0.84)) if finite.size else "",
                        "mean_tree_level_mass_ratio": float(np.mean(finite)) if finite.size else "",
                        "n_zero_mass_selections": sum(v == 0 for v in tree_ratios),
                    })
    write_csv(args.output / "r2_mass_recovery_per_tree.csv", mass_recovery_tree)
    write_csv(args.output / "r2_mass_recovery_summary.csv", mass_recovery_summary)

    reference_growth = []
    for band in ("V", "g"):
        for condition in CONDITIONS:
            for metric, field in (("half_light_radius", "central_bound_R50_pkpc"),
                                  ("half_mass_radius", "central_bound_R50_mass_pkpc")):
                tree_growth = []
                for rank in cohort:
                    los_values = []
                    for los in LOS:
                        start = as_float(get_row(rank, ENDPOINTS[0], band,
                                                 "central_bound_stars_reference", los, condition), field)
                        end = as_float(get_row(rank, ENDPOINTS[1], band,
                                               "central_bound_stars_reference", los, condition), field)
                        if min(start, end) > 0 and np.all(np.isfinite([start, end])):
                            los_values.append(float(np.log10(end / start)))
                    if los_values:
                        tree_growth.append(float(np.median(los_values)))
                stats = interval(tree_growth, rng)
                reference_growth.append({
                    "band": band, "observer_condition": condition, "reference_metric": metric,
                    "n_trees": stats["n_trees_measurable"],
                    "median_reference_growth_dex": stats["median_bias_dex"],
                    "median_ci95": json.dumps(stats["median_bias_ci95"]),
                    "mean_reference_growth_dex": stats["mean_bias_dex"],
                    "mean_ci95": json.dumps(stats["mean_bias_ci95"]),
                    "p16_tree_growth_dex": stats["p16_bias_dex"],
                    "p84_tree_growth_dex": stats["p84_bias_dex"],
                })
    write_csv(args.output / "r3_reference_size_growth.csv", reference_growth)

    metadata = {
        "simulation": "TNG300-1",
        "n_descendant_branches": len(cohort),
        "endpoint_detection_definition": (
            "central image pixel lies inside the connected isophotal region in all three projections "
            "at both snapshots 50 and 84, separately by band and observer condition"),
        "full_cohort_failure_handling": (
            "all 36 branches retained in denominators; zero or undefined endpoint selections are "
            "counted and reported; growth statistics use branches with measurable positive masses/radii"),
        "physical_endpoint_redshifts": {str(k): v for k, v in SNAP_TO_Z.items()},
        "methods": METHODS,
        "bootstrap_replicates": N_BOOT,
        "bootstrap_seed": SEED,
        "score_directory": args.scores.resolve().relative_to(HERE).as_posix(),
        "output_directory": args.output.resolve().relative_to(HERE).as_posix(),
    }
    (args.output / "analysis_metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")
    print(json.dumps({"cohort": len(cohort), "endpoint_detected_trees": {
        f"{band}_{condition}": len(ranks)
        for (band, condition), ranks in endpoint_valid_ranks.items()},
        "output_directory": str(args.output)}, indent=2), flush=True)


if __name__ == "__main__":
    main()
