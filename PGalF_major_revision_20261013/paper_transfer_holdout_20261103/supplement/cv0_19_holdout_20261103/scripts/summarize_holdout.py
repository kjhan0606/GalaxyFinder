#!/usr/bin/env python3
"""Summarize the locked CV0--19 holdout with hosts as statistical units."""
from __future__ import annotations

import argparse
import csv
import json
import math
from collections import defaultdict
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parent
TAG = "cv019_holdout_20261103"
CALIBRATION = {"V": "0.12", "g": "0.10"}
BASELINE = ("FoF", "0.20")
INVALID_REASONS = ("fallback", "empty_contour", "zero_contour_light", "edge_touch")


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream))


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        path.write_text("\n")
        return
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def number(row: dict[str, str], key: str) -> float:
    value = row.get(key, "")
    try:
        return float(value)
    except (TypeError, ValueError):
        return float("nan")


def median(rows: list[dict[str, object]], key: str) -> float:
    values = np.asarray([float(row[key]) for row in rows], dtype=float)
    values = values[np.isfinite(values)]
    return float(np.median(values)) if values.size else float("nan")


def mean_ci(values: np.ndarray, rng: np.random.Generator, draws: int = 2000) -> tuple[float, float]:
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if not values.size:
        return float("nan"), float("nan")
    sample_idx = rng.integers(0, values.size, size=(draws, values.size))
    estimates = np.mean(values[sample_idx], axis=1)
    return tuple(float(x) for x in np.quantile(estimates, [0.025, 0.975]))


def box_host_interval(
    keyed_values: list[tuple[int, float]], statistic: str,
    rng: np.random.Generator, draws: int = 2000,
) -> tuple[float, float]:
    by_box: dict[int, list[float]] = defaultdict(list)
    for box, value in keyed_values:
        if np.isfinite(value):
            by_box[int(box)].append(float(value))
    by_box_array = {box: np.asarray(values) for box, values in by_box.items()}
    boxes = np.asarray(sorted(by_box_array), dtype=int)
    if not boxes.size:
        return float("nan"), float("nan")
    box_values = [by_box_array[int(box)] for box in boxes]
    counts = np.asarray([len(values) for values in box_values], dtype=int)
    max_count = int(np.max(counts))
    padded = np.full((len(boxes), max_count), np.nan, dtype=float)
    for index, values in enumerate(box_values):
        padded[index, :len(values)] = values
    sampled_boxes = rng.integers(0, len(boxes), size=(draws, len(boxes)))
    samples = np.full((draws, len(boxes) * max_count), np.nan, dtype=float)
    within = np.arange(max_count)[None, :]
    for column in range(len(boxes)):
        selected = sampled_boxes[:, column]
        selected_counts = counts[selected]
        indices = (rng.random((draws, max_count)) * selected_counts[:, None]).astype(int)
        values = padded[selected[:, None], indices]
        values[within >= selected_counts[:, None]] = np.nan
        start = column * max_count
        samples[:, start:start + max_count] = values
    if statistic == "mean":
        estimates = np.nanmean(samples, axis=1)
    elif statistic == "p90":
        estimates = np.nanquantile(samples, 0.90, axis=1)
    else:
        raise ValueError(f"Unknown bootstrap statistic: {statistic}")
    return tuple(float(x) for x in np.quantile(estimates, [0.025, 0.975]))


def group_names(row: dict[str, object], band: str) -> list[str]:
    if row["method"] == "aperture30pkpc":
        return ["aperture30pkpc"]
    b = str(row["b_star"])
    names = []
    if b == "0.10":
        names.append("fixed_b010")
    if b == CALIBRATION[band]:
        names.append("darwin_calibrated")
    if b == "0.20":
        names.append("baseline_b020")
    return names


def host_table(projection_rows: list[dict[str, object]]) -> list[dict[str, object]]:
    groups: dict[tuple[int, int, str, str], list[dict[str, object]]] = defaultdict(list)
    for row in projection_rows:
        groups[(int(row["cv_box"]), int(row["fof_group_index"]), str(row["band"]), str(row["method_key"]))].append(row)
    result = []
    for (box, group, band, method), rows in sorted(groups.items()):
        if len(rows) != 3:
            raise ValueError(f"Expected 3 views for CV_{box} group {group} {band} {method}; found {len(rows)}")
        invalid = any(bool(row["invalid_target"]) for row in rows)
        row = {
            "cv_box": box,
            "fof_group_index": group,
            "band": band,
            "method_key": method,
            "b_star": str(rows[0]["b_star"]),
            "n_projections": len(rows),
            "n_invalid_contours": sum(bool(r["invalid_target"]) for r in rows),
            "invalid_contour_any_view": int(invalid),
        }
        metric_columns = (
            "iou", "completeness", "precision", "contamination",
            "log10_luminosity_ratio", "abs_log10_luminosity_error",
            "log10_R50_ratio", "abs_log10_R50_error",
        )
        for metric in metric_columns:
            row[f"median_{metric}"] = median(rows, metric)
        row["fraction_invalid_views"] = row["n_invalid_contours"] / 3.0
        row["iou_below_0p90"] = int(np.isfinite(float(row["median_iou"])) and float(row["median_iou"]) < 0.90)
        log2 = math.log10(2.0)
        row["failure_invalid_contour"] = int(invalid)
        row["failure_contamination"] = int(
            np.isfinite(float(row["median_contamination"]))
            and float(row["median_contamination"]) > 0.20
        )
        row["failure_luminosity"] = int(
            np.isfinite(float(row["median_log10_luminosity_ratio"]))
            and abs(float(row["median_log10_luminosity_ratio"])) > log2
        )
        row["failure_R50"] = int(
            np.isfinite(float(row["median_log10_R50_ratio"]))
            and abs(float(row["median_log10_R50_ratio"])) > log2
        )
        median_failure = any(row[key] for key in (
            "failure_invalid_contour", "failure_contamination",
            "failure_luminosity", "failure_R50",
        ))
        view_failure = any(bool(r["large_failure_projection"]) for r in rows)
        row["large_failure_host_median"] = int(median_failure)
        row["large_failure_any_orientation"] = int(view_failure)
        result.append(row)
    return result


def sample_summary(
    rows: list[dict[str, object]], rng: np.random.Generator,
    include_bootstrap: bool = True,
) -> dict[str, object]:
    out: dict[str, object] = {
        "band": rows[0]["band"],
        "method_key": rows[0]["method_key"],
        "b_star": rows[0]["b_star"],
        "n_hosts": len(rows),
    }
    metrics = (
        "iou", "completeness", "precision", "contamination",
        "log10_luminosity_ratio", "log10_R50_ratio",
    )
    for metric in metrics:
        values = np.asarray([float(r[f"median_{metric}"]) for r in rows], dtype=float)
        values = values[np.isfinite(values)]
        out[f"mean_host_median_{metric}"] = float(np.mean(values)) if values.size else float("nan")
        out[f"median_host_median_{metric}"] = float(np.median(values)) if values.size else float("nan")
        out[f"host_{metric}_p16"] = float(np.quantile(values, 0.16)) if values.size else float("nan")
        out[f"host_{metric}_p84"] = float(np.quantile(values, 0.84)) if values.size else float("nan")
    for metric in ("abs_log10_luminosity_error", "abs_log10_R50_error"):
        values = np.asarray([float(r[f"median_{metric}"]) for r in rows], dtype=float)
        values = values[np.isfinite(values)]
        out[f"n_hosts_valid_{metric}"] = int(values.size)
        out[f"p90_host_{metric}"] = float(np.quantile(values, 0.90)) if values.size else float("nan")
        if include_bootstrap:
            keyed = [(int(r["cv_box"]), float(r[f"median_{metric}"])) for r in rows]
            low, high = box_host_interval(keyed, "p90", rng)
            out[f"p90_host_{metric}_box_host_95_low"] = low
            out[f"p90_host_{metric}_box_host_95_high"] = high
    for label in (
        "iou_below_0p90", "large_failure_host_median",
        "large_failure_any_orientation", "invalid_contour_any_view",
        "failure_invalid_contour", "failure_contamination",
        "failure_luminosity", "failure_R50",
    ):
        values = np.asarray([float(r[label]) for r in rows], dtype=float)
        out[f"fraction_hosts_{label}"] = float(np.mean(values))
        if include_bootstrap:
            keyed = [(int(r["cv_box"]), float(r[label])) for r in rows]
            low, high = box_host_interval(keyed, "mean", rng)
            out[f"fraction_hosts_{label}_box_host_95_low"] = low
            out[f"fraction_hosts_{label}_box_host_95_high"] = high
    for metric in (
        "iou", "completeness", "precision", "contamination",
        "log10_luminosity_ratio", "log10_R50_ratio",
    ):
        if include_bootstrap:
            keyed = [(int(r["cv_box"]), float(r[f"median_{metric}"])) for r in rows]
            low, high = box_host_interval(keyed, "mean", rng)
            out[f"mean_host_median_{metric}_box_host_95_low"] = low
            out[f"mean_host_median_{metric}_box_host_95_high"] = high
    return out


def paired_contrasts(halo_rows: list[dict[str, object]], seed: int = 20261103) -> list[dict[str, object]]:
    by_key = {
        (int(row["cv_box"]), int(row["fof_group_index"]), str(row["band"]), str(row["method_key"])): row
        for row in halo_rows
    }
    rng = np.random.default_rng(seed)
    contrasts = []
    metrics = (
        "iou", "completeness", "precision", "contamination",
        "log10_luminosity_ratio", "abs_log10_luminosity_error",
        "log10_R50_ratio", "abs_log10_R50_error",
        "large_failure_host_median", "large_failure_any_orientation",
        "failure_invalid_contour", "failure_contamination",
        "failure_luminosity", "failure_R50",
    )
    treatments = ("fixed_b010", "darwin_calibrated", "aperture30pkpc")
    for band in ("V", "g"):
        host_ids = sorted({(int(r["cv_box"]), int(r["fof_group_index"])) for r in halo_rows if r["band"] == band})
        for treatment in treatments:
            pairs = []
            for box, group in host_ids:
                base = by_key.get((box, group, band, "baseline_b020"))
                comp = by_key.get((box, group, band, treatment))
                if base is None or comp is None:
                    raise ValueError(f"Missing paired host {box}/{group}/{band}/{treatment}")
                pairs.append((box, base, comp))
            for metric in metrics:
                values = []
                keyed = []
                column = metric if metric.startswith(("large_failure_", "failure_")) else f"median_{metric}"
                for box, base, comp in pairs:
                    delta = float(comp[column]) - float(base[column])
                    if np.isfinite(delta):
                        values.append(delta)
                        keyed.append((box, delta))
                low, high = mean_ci(np.asarray(values), rng)
                box_low, box_high = box_host_interval(keyed, "mean", rng)
                contrasts.append({
                    "band": band,
                    "baseline_method": "baseline_b020",
                    "comparison_method": treatment,
                    "metric": metric,
                    "n_paired_hosts": len(values),
                    "mean_paired_host_difference": float(np.mean(values)) if values else float("nan"),
                    "host_bootstrap_95_low": low,
                    "host_bootstrap_95_high": high,
                    "box_then_host_bootstrap_95_low": box_low,
                    "box_then_host_bootstrap_95_high": box_high,
                })
        # Compare each FoF rule directly with the same-host 30 pkpc aperture.
        for treatment in ("fixed_b010", "darwin_calibrated"):
            pairs = []
            for box, group in host_ids:
                aperture = by_key.get((box, group, band, "aperture30pkpc"))
                comp = by_key.get((box, group, band, treatment))
                if aperture is None or comp is None:
                    raise ValueError(f"Missing aperture-paired host {box}/{group}/{band}/{treatment}")
                pairs.append((box, aperture, comp))
            for metric in metrics:
                values, keyed = [], []
                column = metric if metric.startswith(("large_failure_", "failure_")) else f"median_{metric}"
                for box, aperture, comp in pairs:
                    delta = float(comp[column]) - float(aperture[column])
                    if np.isfinite(delta):
                        values.append(delta)
                        keyed.append((box, delta))
                low, high = mean_ci(np.asarray(values), rng)
                box_low, box_high = box_host_interval(keyed, "mean", rng)
                contrasts.append({
                    "band": band,
                    "baseline_method": "aperture30pkpc",
                    "comparison_method": treatment,
                    "metric": metric,
                    "n_paired_hosts": len(values),
                    "mean_paired_host_difference": float(np.mean(values)) if values else float("nan"),
                    "host_bootstrap_95_low": low,
                    "host_bootstrap_95_high": high,
                    "box_then_host_bootstrap_95_low": box_low,
                    "box_then_host_bootstrap_95_high": box_high,
                })
        # The V-band Darwin value is 0.12, while the common fixed comparison is 0.10.
        # In g both definitions are the same method and no duplicate contrast is emitted.
        if band == "V":
            pairs = []
            for box, group in host_ids:
                fixed = by_key[(box, group, band, "fixed_b010")]
                calibrated = by_key[(box, group, band, "darwin_calibrated")]
                pairs.append((box, fixed, calibrated))
            for metric in metrics:
                values, keyed = [], []
                column = metric if metric.startswith(("large_failure_", "failure_")) else f"median_{metric}"
                for box, fixed, calibrated in pairs:
                    delta = float(calibrated[column]) - float(fixed[column])
                    if np.isfinite(delta):
                        values.append(delta)
                        keyed.append((box, delta))
                low, high = mean_ci(np.asarray(values), rng)
                box_low, box_high = box_host_interval(keyed, "mean", rng)
                contrasts.append({
                    "band": band,
                    "baseline_method": "fixed_b010",
                    "comparison_method": "darwin_calibrated",
                    "metric": metric,
                    "n_paired_hosts": len(values),
                    "mean_paired_host_difference": float(np.mean(values)) if values else float("nan"),
                    "host_bootstrap_95_low": low,
                    "host_bootstrap_95_high": high,
                    "box_then_host_bootstrap_95_low": box_low,
                    "box_then_host_bootstrap_95_high": box_high,
                })
    return contrasts


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", type=Path, default=ROOT / "sample_manifest.csv")
    parser.add_argument("--results", type=Path, default=ROOT / "results")
    parser.add_argument("--tag", default=TAG)
    parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args()
    if args.outdir.exists() and any(args.outdir.iterdir()):
        raise SystemExit(f"Refusing to overwrite nonempty output directory: {args.outdir}")

    manifest = read_csv(args.manifest)
    if len(manifest) != 130:
        raise ValueError(f"Expected the frozen 130-host manifest, got {len(manifest)}")
    projection_rows: list[dict[str, object]] = []
    missing = []
    for host in manifest:
        box = int(host["cv_box"])
        group = int(host["fof_group_index"])
        path = args.results / f"CV{box:02d}_group{group:05d}_{args.tag}" / "metrics.csv"
        if not path.is_file():
            missing.append({"cv_box": box, "fof_group_index": group, "missing": str(path)})
            continue
        metrics = read_csv(path)
        targets = {
            (row["line_of_sight"], row["band"]): row
            for row in metrics
            if row["target_variant"] == "potential" and row["method"] == "aperture30pkpc"
        }
        for raw in metrics:
            if raw["target_variant"] != "potential":
                continue
            target = targets.get((raw["line_of_sight"], raw["band"]))
            if target is None:
                raise ValueError(f"Missing target row for CV_{box} group {group}: {raw['line_of_sight']} {raw['band']}")
            fallback = number(target, "target_component_fallback") != 0
            empty = number(target, "target_npixel") == 0
            zero_light = number(target, "L_contour") <= 0
            edge = number(target, "target_touches_image_edge") == 1
            reasons = []
            if fallback:
                reasons.append("fallback")
            if empty:
                reasons.append("empty_contour")
            if zero_light:
                reasons.append("zero_contour_light")
            if edge:
                reasons.append("edge_touch")
            invalid = bool(reasons)
            method = raw["method"]
            b_star = "aperture" if method == "aperture30pkpc" else f"{number(raw, 'b_star'):.2f}"
            values = {}
            for key in ("iou", "completeness", "precision"):
                values[key] = number(raw, key)
            log_lum = log_r50 = float("nan")
            selected_lum = number(raw, "L_selected")
            contour_lum = number(raw, "L_contour")
            r50_selected = number(raw, "R50_selected_kpc")
            r50_contour = number(raw, "R50_contour_kpc")
            if not invalid and selected_lum > 0 and contour_lum > 0:
                log_lum = math.log10(selected_lum / contour_lum)
            if not invalid and r50_selected > 0 and r50_contour > 0:
                log_r50 = math.log10(r50_selected / r50_contour)
            if invalid:
                values["iou"] = 0.0
                values["completeness"] = 0.0
                values["precision"] = float("nan")
            contamination = 1.0 - values["precision"] if np.isfinite(values["precision"]) else float("nan")
            projection_large_failure = (
                invalid
                or (np.isfinite(log_lum) and abs(log_lum) > math.log10(2.0))
                or (np.isfinite(log_r50) and abs(log_r50) > math.log10(2.0))
                or (np.isfinite(contamination) and contamination > 0.20)
            )
            method_keys = group_names({"method": method, "b_star": b_star}, raw["band"])
            for method_key in method_keys:
                projection_rows.append({
                    "cv_box": box,
                    "fof_group_index": group,
                    "band": raw["band"],
                    "line_of_sight": raw["line_of_sight"],
                    "method": method,
                    "method_key": method_key,
                    "b_star": b_star,
                    "invalid_target": int(invalid),
                    "invalid_reasons": ";".join(reasons),
                    "iou": values["iou"],
                    "completeness": values["completeness"],
                    "precision": values["precision"],
                    "contamination": contamination,
                    "log10_luminosity_ratio": log_lum,
                    "abs_log10_luminosity_error": abs(log_lum) if np.isfinite(log_lum) else float("nan"),
                    "log10_R50_ratio": log_r50,
                    "abs_log10_R50_error": abs(log_r50) if np.isfinite(log_r50) else float("nan"),
                    "large_failure_projection": int(projection_large_failure),
                })
    if missing:
        raise FileNotFoundError(f"Missing measurements for {len(missing)} frozen hosts; first: {missing[:3]}")
    halo_rows = host_table(projection_rows)
    summary_rng = np.random.default_rng(20261103)
    summary_rows = []
    groups: dict[tuple[str, str], list[dict[str, object]]] = defaultdict(list)
    for row in halo_rows:
        groups[(str(row["band"]), str(row["method_key"]))].append(row)
    for key in sorted(groups):
        summary_rows.append(sample_summary(groups[key], summary_rng))

    box_groups: dict[tuple[int, str, str], list[dict[str, object]]] = defaultdict(list)
    for row in halo_rows:
        box_groups[(int(row["cv_box"]), str(row["band"]), str(row["method_key"]))].append(row)
    box_rows = []
    for (box, band, method), rows in sorted(box_groups.items()):
        summary = sample_summary(
            rows, np.random.default_rng(20261103 + box), include_bootstrap=False
        )
        box_rows.append({"cv_box": box, **summary})

    contrast_rows = paired_contrasts(halo_rows)
    args.outdir.mkdir(parents=True, exist_ok=True)
    write_csv(args.outdir / "projection_metrics.csv", projection_rows)
    write_csv(args.outdir / "host_metrics.csv", halo_rows)
    write_csv(args.outdir / "box_summary.csv", box_rows)
    write_csv(args.outdir / "paired_contrasts.csv", contrast_rows)
    write_csv(args.outdir / "summary.csv", summary_rows)
    report = {
        "expected_hosts": len(manifest),
        "hosts_with_metrics": len({(r["cv_box"], r["fof_group_index"]) for r in projection_rows}),
        "missing_hosts": missing,
        "target_variant": "potential",
        "calibration_b_star": CALIBRATION,
        "baseline": {"method": "FoF", "b_star": 0.20},
        "comparisons": ["fixed b_star=0.10", "Darwin-calibrated b_star by band", "30 pkpc aperture"],
        "statistical_unit": "host halo; the three projections are summarized within each host",
        "invalid_contour": "fallback from the potential-centred component, no contour pixels, zero contour luminosity, or contour touching the image edge",
        "large_failure_host_median": "any invalid projected contour, median contamination >0.20, or absolute host-median log10 luminosity or R50 ratio > log10(2)",
        "large_failure_any_orientation": "at least one projection has any invalid contour, contamination >0.20, or an absolute log10 luminosity or R50 ratio > log10(2)",
        "absolute_error_quantiles": "the 90th percentile is conditional on valid target contours; invalid contours are reported in both large-failure rates",
        "bootstrap": "2,000 percentile draws; resample CV boxes first, then hosts within each sampled box; paired contrasts use host-matched method differences",
        "summary_rows": summary_rows,
    }
    (args.outdir / "summary.json").write_text(json.dumps(report, indent=2, allow_nan=True) + "\n")
    print(json.dumps({
        "hosts": report["hosts_with_metrics"],
        "projection_rows": len(projection_rows),
        "host_method_rows": len(halo_rows),
        "summary_rows": len(summary_rows),
        "paired_contrasts": len(contrast_rows),
        "output_directory": str(args.outdir),
    }, indent=2))


if __name__ == "__main__":
    main()
