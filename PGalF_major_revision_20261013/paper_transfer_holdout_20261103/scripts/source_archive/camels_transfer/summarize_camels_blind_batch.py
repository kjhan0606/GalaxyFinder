#!/usr/bin/env python3
"""Summarize the 49-halo blind-transfer metrics by host and photometric band."""
from __future__ import annotations

import argparse
import csv
import json
import math
from collections import defaultdict
from pathlib import Path

import numpy as np


REVIEW = Path(__file__).resolve().parent
MANIFEST = REVIEW / "sample_manifest.csv"
RESULTS = REVIEW / "results"
TAG = "blind20261029_v2"
TARGET_VARIANT = "light_peak"


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


def finite_values(rows: list[dict[str, object]], key: str) -> np.ndarray:
    values = np.asarray([float(r[key]) for r in rows], dtype=float)
    return values[np.isfinite(values)]


def host_bootstrap_interval(values: np.ndarray, statistic: str, rng: np.random.Generator, draws: int = 2000) -> tuple[float, float]:
    values = values[np.isfinite(values)]
    if not len(values):
        return float("nan"), float("nan")
    indices = rng.integers(0, len(values), size=(draws, len(values)))
    sample = values[indices]
    if statistic == "mean":
        estimates = np.mean(sample, axis=1)
    elif statistic == "p90":
        estimates = np.quantile(sample, 0.90, axis=1)
    else:
        raise ValueError(f"Unknown bootstrap statistic {statistic}")
    low, high = np.quantile(estimates, [0.025, 0.975])
    return float(low), float(high)


def box_host_bootstrap_interval(
    keyed_values: list[tuple[int, float]], statistic: str,
    rng: np.random.Generator, draws: int = 2000,
) -> tuple[float, float]:
    """Resample CV boxes first, then hosts within each sampled box."""
    by_box: dict[int, np.ndarray] = defaultdict(list)
    for box, value in keyed_values:
        if np.isfinite(value):
            by_box[int(box)].append(float(value))
    by_box = {box: np.asarray(values, dtype=float) for box, values in by_box.items()}
    boxes = np.asarray(sorted(by_box), dtype=int)
    if boxes.size == 0:
        return float("nan"), float("nan")
    estimates = np.empty(draws, dtype=float)
    for draw in range(draws):
        sampled_boxes = rng.choice(boxes, size=boxes.size, replace=True)
        values = []
        for box in sampled_boxes:
            hosts = by_box[int(box)]
            values.extend(rng.choice(hosts, size=hosts.size, replace=True))
        sample = np.asarray(values, dtype=float)
        estimates[draw] = np.mean(sample) if statistic == "mean" else np.quantile(sample, 0.90)
    return tuple(float(value) for value in np.quantile(estimates, [0.025, 0.975]))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", type=Path, default=MANIFEST)
    parser.add_argument("--results", type=Path, default=RESULTS)
    parser.add_argument("--tag", default=TAG)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--target-variant", choices=("light_peak", "center", "potential"), default=TARGET_VARIANT)
    args = parser.parse_args()
    if args.outdir.exists() and any(args.outdir.iterdir()):
        raise SystemExit(f"Refusing to overwrite nonempty output directory: {args.outdir}")

    manifest = read_csv(args.manifest)
    projection_rows = []
    missing = []
    for halo in manifest:
        box = int(halo["cv_box"])
        group = int(halo["fof_group_index"])
        folder = args.results / f"CV{box:02d}_group{group:05d}_{args.tag}"
        metrics_path = folder / "metrics.csv"
        if not metrics_path.is_file():
            missing.append({"cv_box": box, "fof_group_index": group, "missing": str(metrics_path)})
            continue
        metrics = read_csv(metrics_path)
        target_rows = {
            (r["line_of_sight"], r["band"]): r
            for r in metrics
            if r["target_variant"] == args.target_variant and r["method"] == "aperture30pkpc"
        }
        for row in metrics:
            if row["target_variant"] != args.target_variant:
                continue
            key = (row["line_of_sight"], row["band"])
            target = target_rows.get(key)
            if target is None:
                raise ValueError(f"Missing target row for CV{box} group {group} {key}")
            failure = (
                int(float(target["target_component_fallback"])) != 0
                or int(float(target["target_npixel"])) == 0
                or float(target["L_contour"]) <= 0.0
            )
            method = row["method"]
            b_value = "aperture30pkpc" if method == "aperture30pkpc" else f"{float(row['b_star']):.2f}"
            lum_ratio = (
                1.0 + float(row["luminosity_ratio_minus_one"])
                if row["luminosity_ratio_minus_one"] != "" else float("nan")
            )
            r50_selected = float(row["R50_selected_kpc"]) if row["R50_selected_kpc"] != "" else float("nan")
            r50_contour = float(row["R50_contour_kpc"]) if row["R50_contour_kpc"] != "" else float("nan")
            log_lum_ratio = math.log10(lum_ratio) if np.isfinite(lum_ratio) and lum_ratio > 0 else float("nan")
            log_r50_ratio = (
                math.log10(r50_selected / r50_contour)
                if np.isfinite(r50_selected) and np.isfinite(r50_contour) and r50_selected > 0 and r50_contour > 0
                else float("nan")
            )
            precision_value = float(row["precision"]) if row["precision"] != "" and np.isfinite(float(row["precision"])) else float("nan")
            item = {
                "cv_box": box,
                "fof_group_index": group,
                "fof_group_id_raw": int(halo["fof_group_id_raw"]),
                "line_of_sight": row["line_of_sight"],
                "band": row["band"],
                "target_variant": args.target_variant,
                "method": method,
                "b_star": b_value,
                "target_failed": int(failure),
                "target_fallback": int(float(target["target_component_fallback"]) != 0),
                "target_npixel": int(float(target["target_npixel"])),
                "iou": float(row["iou"]) if row["iou"] != "" and np.isfinite(float(row["iou"])) else float("nan"),
                "completeness": float(row["completeness"]) if row["completeness"] != "" and np.isfinite(float(row["completeness"])) else float("nan"),
                "precision": precision_value,
                "contamination_fraction_outside_contour": 1.0 - precision_value if np.isfinite(precision_value) else float("nan"),
                "luminosity_ratio_minus_one": float(row["luminosity_ratio_minus_one"]) if row["luminosity_ratio_minus_one"] != "" and np.isfinite(float(row["luminosity_ratio_minus_one"])) else float("nan"),
                "log10_luminosity_ratio": log_lum_ratio,
                "abs_log10_luminosity_error": abs(log_lum_ratio) if np.isfinite(log_lum_ratio) else float("nan"),
                "R50_selected_kpc": r50_selected,
                "R50_contour_kpc": r50_contour,
                "log10_R50_ratio": log_r50_ratio,
                "abs_log10_R50_error": abs(log_r50_ratio) if np.isfinite(log_r50_ratio) else float("nan"),
            }
            if failure:
                item["iou"] = 0.0
                item["completeness"] = 0.0
                # A missing or fallback contour has no defined luminosity or
                # size error. Keep it in the explicit failure-rate count.
                for key in (
                    "precision", "contamination_fraction_outside_contour",
                    "luminosity_ratio_minus_one", "log10_luminosity_ratio",
                    "abs_log10_luminosity_error", "R50_contour_kpc",
                    "log10_R50_ratio", "abs_log10_R50_error",
                ):
                    item[key] = float("nan")
            projection_rows.append(item)

    # The projection is a repeated view, not an independent statistical unit.
    grouped: dict[tuple, list[dict[str, object]]] = defaultdict(list)
    for row in projection_rows:
        key = (row["cv_box"], row["fof_group_index"], row["band"], row["method"], row["b_star"])
        grouped[key].append(row)
    halo_rows = []
    for (box, group, band, method, b_value), rows in grouped.items():
        failures = sum(int(r["target_failed"]) for r in rows)
        row = {
            "cv_box": box, "fof_group_index": group, "band": band,
            "method": method, "b_star": b_value,
            "n_projections": len(rows), "n_failed_projections": failures,
            "fraction_failed_projections": failures / len(rows) if rows else float("nan"),
        }
        for metric in ("iou", "completeness", "precision", "contamination_fraction_outside_contour", "luminosity_ratio_minus_one", "log10_luminosity_ratio", "abs_log10_luminosity_error", "R50_selected_kpc", "R50_contour_kpc", "log10_R50_ratio", "abs_log10_R50_error"):
            values = finite_values(rows, metric)
            row[f"median_{metric}"] = float(np.median(values)) if len(values) else float("nan")
        halo_rows.append(row)

    box_groups: dict[tuple, list[dict[str, object]]] = defaultdict(list)
    for row in halo_rows:
        box_groups[(row["cv_box"], row["band"], row["method"], row["b_star"])].append(row)
    box_summary_rows = []
    for (box, band, method, b_value), rows in box_groups.items():
        box_lum_bias = finite_values(rows, "median_log10_luminosity_ratio")
        box_lum_error = finite_values(rows, "median_abs_log10_luminosity_error")
        box_r50_bias = finite_values(rows, "median_log10_R50_ratio")
        box_r50_error = finite_values(rows, "median_abs_log10_R50_error")
        item = {
            "cv_box": box, "band": band, "method": method, "b_star": b_value,
            "n_halos": len(rows),
            "mean_halo_median_iou": float(np.mean([float(r["median_iou"]) for r in rows])),
            "median_halo_median_iou": float(np.median([float(r["median_iou"]) for r in rows])),
            "fraction_halos_iou_below_0p90": float(np.mean([float(r["median_iou"]) < 0.90 for r in rows])),
            "mean_halo_contamination_fraction_outside_contour": float(np.mean(finite_values(rows, "median_contamination_fraction_outside_contour"))) if len(finite_values(rows, "median_contamination_fraction_outside_contour")) else float("nan"),
            "mean_halo_signed_log10_luminosity_bias": float(np.mean(box_lum_bias)) if len(box_lum_bias) else float("nan"),
            "p90_halo_abs_log10_luminosity_error": float(np.quantile(box_lum_error, 0.90)) if len(box_lum_error) else float("nan"),
            "mean_halo_signed_log10_R50_bias": float(np.mean(box_r50_bias)) if len(box_r50_bias) else float("nan"),
            "p90_halo_abs_log10_R50_error": float(np.quantile(box_r50_error, 0.90)) if len(box_r50_error) else float("nan"),
        }
        box_summary_rows.append(item)

    # Paired host contrasts preserve the same simulation hosts across methods.
    contrast_rows = []
    contrast_rng = np.random.default_rng(20261030)
    baseline_key = ("FoF", "0.20")
    treatments = (("FoF", "0.10"), ("FoF", "0.09"), ("aperture30pkpc", "aperture30pkpc"))
    contrast_metrics = {
        "iou": "median_iou",
        "completeness": "median_completeness",
        "contamination_fraction_outside_contour": "median_contamination_fraction_outside_contour",
        "abs_log10_luminosity_error": "median_abs_log10_luminosity_error",
        "abs_log10_R50_error": "median_abs_log10_R50_error",
    }
    for band in sorted({str(row["band"]) for row in halo_rows}):
        baseline_hosts = {
            (int(row["cv_box"]), int(row["fof_group_index"])): row
            for row in halo_rows
            if str(row["band"]) == band and (str(row["method"]), str(row["b_star"])) == baseline_key
        }
        for treatment_method, treatment_b in treatments:
            treatment_hosts = {
                (int(row["cv_box"]), int(row["fof_group_index"])): row
                for row in halo_rows
                if str(row["band"]) == band
                and (str(row["method"]), str(row["b_star"])) == (treatment_method, treatment_b)
            }
            common_hosts = sorted(set(baseline_hosts) & set(treatment_hosts))
            if len(common_hosts) != len(manifest):
                raise ValueError(f"Paired method comparison has {len(common_hosts)}/{len(manifest)} hosts for {band} {treatment_method} {treatment_b}")
            deltas = {}
            for label, column in contrast_metrics.items():
                deltas[label] = np.asarray([
                    float(treatment_hosts[key][column]) - float(baseline_hosts[key][column])
                    for key in common_hosts
                ], dtype=float)
            deltas["fraction_halos_iou_below_0p90"] = np.asarray([
                float(float(treatment_hosts[key]["median_iou"]) < 0.90)
                - float(float(baseline_hosts[key]["median_iou"]) < 0.90)
                for key in common_hosts
            ], dtype=float)
            box_keyed_hosts = [
                (key, baseline_hosts[key], treatment_hosts[key]) for key in common_hosts
            ]
            for metric, values in deltas.items():
                values = values[np.isfinite(values)]
                low, high = host_bootstrap_interval(values, "mean", contrast_rng)
                metric_column = contrast_metrics.get(metric)
                if metric_column:
                    box_keyed = [
                        (int(base_row["cv_box"]), float(treat_row[metric_column]) - float(base_row[metric_column]))
                        for _key, base_row, treat_row in box_keyed_hosts
                        if np.isfinite(float(treat_row[metric_column]))
                        and np.isfinite(float(base_row[metric_column]))
                    ]
                else:
                    box_keyed = [
                        (int(base_row["cv_box"]),
                         float(float(treat_row["median_iou"]) < 0.90)
                         - float(float(base_row["median_iou"]) < 0.90))
                        for _key, base_row, treat_row in box_keyed_hosts
                    ]
                box_low, box_high = box_host_bootstrap_interval(box_keyed, "mean", contrast_rng)
                contrast_rows.append({
                    "band": band,
                    "baseline_method": baseline_key[0], "baseline_b_star": baseline_key[1],
                    "comparison_method": treatment_method, "comparison_b_star": treatment_b,
                    "metric": metric, "n_paired_hosts": len(values),
                    "mean_paired_host_difference": float(np.mean(values)) if len(values) else float("nan"),
                    "paired_host_bootstrap_95_low": low,
                    "paired_host_bootstrap_95_high": high,
                    "paired_box_then_host_bootstrap_95_low": box_low,
                    "paired_box_then_host_bootstrap_95_high": box_high,
                })

    summary_rows = []
    summary_groups: dict[tuple, list[dict[str, object]]] = defaultdict(list)
    bootstrap_rng = np.random.default_rng(20261029)
    for row in halo_rows:
        summary_groups[(row["band"], row["method"], row["b_star"])].append(row)
    for (band, method, b_value), rows in summary_groups.items():
        summary = {"band": band, "method": method, "b_star": b_value, "n_halos": len(rows)}
        summary["median_failure_fraction"] = float(np.median([r["fraction_failed_projections"] for r in rows]))
        for metric in ("iou", "completeness", "precision", "contamination_fraction_outside_contour", "luminosity_ratio_minus_one", "log10_luminosity_ratio", "abs_log10_luminosity_error", "R50_selected_kpc", "R50_contour_kpc", "log10_R50_ratio", "abs_log10_R50_error"):
            values = finite_values(rows, f"median_{metric}")
            summary[f"median_host_{metric}"] = float(np.median(values)) if len(values) else float("nan")
            if metric == "iou":
                summary["mean_host_median_iou"] = float(np.mean(values)) if len(values) else float("nan")
            summary[f"host_{metric}_p16"] = float(np.quantile(values, 0.16)) if len(values) else float("nan")
            summary[f"host_{metric}_p84"] = float(np.quantile(values, 0.84)) if len(values) else float("nan")
            if metric in ("abs_log10_luminosity_error", "abs_log10_R50_error"):
                summary[f"host_{metric}_p90"] = float(np.quantile(values, 0.90)) if len(values) else float("nan")
        signed_lum = finite_values(rows, "median_log10_luminosity_ratio")
        signed_r50 = finite_values(rows, "median_log10_R50_ratio")
        summary["mean_host_signed_log10_luminosity_bias"] = float(np.mean(signed_lum)) if len(signed_lum) else float("nan")
        summary["mean_host_signed_log10_R50_bias"] = float(np.mean(signed_r50)) if len(signed_r50) else float("nan")
        summary["n_halos_iou_below_0p90"] = sum(float(r["median_iou"]) < 0.90 for r in rows if np.isfinite(float(r["median_iou"])))
        summary["fraction_halos_iou_below_0p90"] = summary["n_halos_iou_below_0p90"] / len(rows) if rows else float("nan")
        summary["n_halos_with_any_failed_projection"] = sum(float(r["fraction_failed_projections"]) > 0 for r in rows)
        summary["fraction_halos_with_any_failed_projection"] = summary["n_halos_with_any_failed_projection"] / len(rows) if rows else float("nan")
        bootstrap_targets = {
            "mean_host_median_iou": (finite_values(rows, "median_iou"), "mean"),
            "mean_host_signed_log10_luminosity_bias": (signed_lum, "mean"),
            "host_abs_log10_luminosity_error_p90": (finite_values(rows, "median_abs_log10_luminosity_error"), "p90"),
            "host_abs_log10_R50_error_p90": (finite_values(rows, "median_abs_log10_R50_error"), "p90"),
            "fraction_halos_iou_below_0p90": (
                np.asarray([float(r["median_iou"]) < 0.90 for r in rows], dtype=float), "mean"
            ),
        }
        for metric_name, (values, statistic) in bootstrap_targets.items():
            low, high = host_bootstrap_interval(values, statistic, bootstrap_rng)
            summary[f"{metric_name}_host_bootstrap_95_low"] = low
            summary[f"{metric_name}_host_bootstrap_95_high"] = high
            column = {
                "mean_host_median_iou": "median_iou",
                "mean_host_signed_log10_luminosity_bias": "median_log10_luminosity_ratio",
                "host_abs_log10_luminosity_error_p90": "median_abs_log10_luminosity_error",
                "host_abs_log10_R50_error_p90": "median_abs_log10_R50_error",
                "fraction_halos_iou_below_0p90": "median_iou",
            }[metric_name]
            keyed = []
            for host in rows:
                value = float(host[column])
                if metric_name == "fraction_halos_iou_below_0p90":
                    value = float(value < 0.90)
                keyed.append((int(host["cv_box"]), value))
            box_low, box_high = box_host_bootstrap_interval(keyed, statistic, bootstrap_rng)
            summary[f"{metric_name}_box_then_host_bootstrap_95_low"] = box_low
            summary[f"{metric_name}_box_then_host_bootstrap_95_high"] = box_high
        summary_rows.append(summary)

    args.outdir.mkdir(parents=True, exist_ok=True)
    write_csv(args.outdir / "projection_metrics.csv", projection_rows)
    write_csv(args.outdir / "halo_metrics.csv", halo_rows)
    write_csv(args.outdir / "box_summary.csv", box_summary_rows)
    write_csv(args.outdir / "paired_host_contrasts.csv", contrast_rows)
    write_csv(args.outdir / "summary.csv", summary_rows)
    report = {
        "expected_halos": len(manifest),
        "halos_with_metrics": len({(r["cv_box"], r["fof_group_index"]) for r in projection_rows}),
        "missing_halos": missing,
        "target_variant": args.target_variant,
        "failed_projection_definition": "fallback from the central pixel, no contour pixels, or zero contour luminosity; failed IoU and completeness are set to zero and failure rates are reported.",
        "large_failure_definition": "halo-level median IoU below 0.90; projections are first summarized within halo.",
        "contamination_definition": "one minus precision: fraction of method-selected stellar luminosity outside the observational contour, weighted by stellar light.",
        "luminosity_bias_and_error": "per-halo median log10(L_selected/L_contour); report equal-host mean signed bias and the 90th percentile of absolute per-halo error.",
        "size_bias_and_error": "per-halo median log10(R50_selected/R50_contour); report equal-host mean signed bias and the 90th percentile of absolute per-halo error.",
        "uncertainty": "95% percentile intervals from 2,000 nested bootstrap draws: resample CV boxes, then hosts within each sampled box; projections are summarized within each host first.",
        "absolute_error_quantiles": "conditional on a valid observational contour; no-contour/fallback cases are counted separately as failures.",
        "statistical_unit": "halo; three projections are summarized within each halo before the across-halo summary.",
        "box_summary": "per-CV-box host summaries are reported separately to expose box-to-box variation; boxes are not pooled as extra halos.",
        "paired_method_contrasts": "b=0.10, b=0.09, and 30 pkpc aperture are contrasted with b=0.20 on the same host set; intervals bootstrap paired host differences.",
        "summary_rows": summary_rows,
    }
    (args.outdir / "summary.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({
        "expected_halos": len(manifest),
        "halos_with_metrics": report["halos_with_metrics"],
        "missing_halos": len(missing),
        "summary_rows": len(summary_rows),
        "box_summary_rows": len(box_summary_rows),
        "paired_contrast_rows": len(contrast_rows),
        "output_directory": str(args.outdir),
    }, indent=2))


if __name__ == "__main__":
    main()
