#!/usr/bin/env python3
"""Leave-one-host-out validation in the consistent z=0.1904 Darwin output."""

from __future__ import annotations

import csv
import hashlib
import json
import os
import platform
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
PAPER = HERE.parents[1]
SUPPLEMENT = PAPER / "supplement"
METADATA = SUPPLEMENT / "halo_metadata.tsv"
CURVES = SUPPLEMENT / "recomputed_curves.tsv"
TARGET_STEP = "01001"
EXPECTED_HOSTS = 12
VIEWS = ("x", "y", "z")
BANDS = ("V", "g")
FIXED_B = (0.10, 0.20)


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open() as stream:
        return list(
            csv.DictReader(
                (line for line in stream if line.strip() and not line.startswith("#")),
                delimiter="\t",
            )
        )


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        raise ValueError(f"Refusing to write an empty table: {path}")
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def median(values: list[float]) -> float:
    return float(np.median(np.asarray(values, dtype=float)))


def main() -> None:
    metadata = read_tsv(METADATA)
    target_rows = [
        row for row in metadata
        if row["step"] == TARGET_STEP
        and int(row["interior_V"])
        and int(row["interior_g"])
    ]
    entries = sorted((row["step"], int(row["nstar"])) for row in target_rows)
    if len(entries) != EXPECTED_HOSTS:
        raise ValueError(f"Expected {EXPECTED_HOSTS} retained hosts, found {len(entries)}")
    z_values = {float(row["z"]) for row in target_rows}
    if len(z_values) != 1:
        raise ValueError(f"Target output has inconsistent redshifts: {z_values}")

    curves: dict[tuple[str, int, str, str], dict[float, tuple[float, float, float]]] = defaultdict(dict)
    for row in read_tsv(CURVES):
        key = (row["step"], int(row["nstar"]))
        if row["variant"] == "consistent_distance" and key in entries:
            curves[(*key, row["band"], row["los"])][float(row["b"])] = (
                float(row["iou"]),
                float(row["completeness"]),
                float(row["precision"]),
            )

    b_grid = np.asarray(
        sorted({b for key in curves for b in curves[key]}), dtype=float
    )
    if not set(FIXED_B).issubset(set(b_grid.tolist())):
        raise ValueError(f"Fixed comparison values are absent from the grid: {b_grid}")

    view_rows: list[dict[str, object]] = []
    host_rows: list[dict[str, object]] = []
    fold_rows: list[dict[str, object]] = []
    summary: dict[str, object] = {
        "analysis": "leave-one-host-out validation within one Darwin output",
        "inputs": {"halo_metadata.tsv": sha256(METADATA), "recomputed_curves.tsv": sha256(CURVES)},
        "sample": {
            "step": TARGET_STEP,
            "redshift": next(iter(z_values)),
            "n_retained_hosts": len(entries),
            "retention": "interior maxima in all three views in both bands",
            "photometry": "consistent_distance, Darwin cosmology",
        },
        "method": {
            "training": "hold out one host; derive each training host's parameter as the median of its three view-specific IoU maxima; take the median across the 11 training hosts",
            "evaluation": "report each held-out projection; form one host value as the median over three projections",
            "fixed_comparisons": list(FIXED_B),
            "tie_rule": "select the smallest tested parameter at a tied maximum",
            "uncertainty": "descriptive only; leave-one-host-out folds overlap and share one output and simulation volume",
        },
        "bands": {},
    }

    for band in BANDS:
        for test_key in entries:
            training_keys = [key for key in entries if key != test_key]
            train_entry_best = []
            for train_key in training_keys:
                view_best = []
                for los in VIEWS:
                    scores = np.asarray(
                        [curves[(*train_key, band, los)][float(b)][0] for b in b_grid]
                    )
                    view_best.append(float(b_grid[int(np.argmax(scores))]))
                train_entry_best.append(median(view_best))
            selected_b = median(train_entry_best)
            if selected_b not in set(b_grid.tolist()):
                raise ValueError(f"Selected parameter {selected_b} is absent from the grid")

            fold_rows.append(
                {
                    "band": band,
                    "held_out_step": test_key[0],
                    "held_out_nstar": test_key[1],
                    "n_training_hosts": len(training_keys),
                    "selected_b": selected_b,
                }
            )
            for los in VIEWS:
                per_b = curves[(*test_key, band, los)]
                row: dict[str, object] = {
                    "band": band,
                    "step": test_key[0],
                    "nstar": test_key[1],
                    "los": los,
                    "n_training_hosts": len(training_keys),
                    "selected_b": selected_b,
                }
                for label, b in (("selected", selected_b), ("fixed_b010", 0.10), ("fixed_b020", 0.20)):
                    values = per_b[b]
                    row[f"{label}_iou"] = values[0]
                    row[f"{label}_completeness"] = values[1]
                    row[f"{label}_precision"] = values[2]
                view_rows.append(row)

        band_views = [row for row in view_rows if row["band"] == band]
        by_host: dict[int, list[dict[str, object]]] = defaultdict(list)
        for row in band_views:
            by_host[int(row["nstar"])].append(row)
        host_metrics: dict[str, list[float]] = defaultdict(list)
        for host_id, rows in sorted(by_host.items()):
            host_row: dict[str, object] = {"band": band, "nstar": host_id}
            for label in ("selected", "fixed_b010", "fixed_b020"):
                for metric in ("iou", "completeness", "precision"):
                    value = median([float(row[f"{label}_{metric}"]) for row in rows])
                    host_row[f"{label}_{metric}"] = value
                    host_metrics[f"{label}_{metric}"].append(value)
            host_row["selected_minus_b020_iou"] = (
                float(host_row["selected_iou"]) - float(host_row["fixed_b020_iou"])
            )
            host_rows.append(host_row)

        gains = [float(row["selected_minus_b020_iou"]) for row in host_rows if row["band"] == band]
        selected_host_iou = host_metrics["selected_iou"]
        selected_min_view_iou = [
            min(float(row["selected_iou"]) for row in by_host[host_id])
            for host_id in by_host
        ]
        selected_b_counts = Counter(float(row["selected_b"]) for row in fold_rows if row["band"] == band)
        band_summary: dict[str, object] = {
            "n_hosts": len(by_host),
            "selected_b_by_fold": {str(b): count for b, count in sorted(selected_b_counts.items())},
            "selected_host_median_iou_median": median(selected_host_iou),
            "selected_host_median_iou_mean": float(np.mean(selected_host_iou)),
            "b020_host_median_iou_median": median(host_metrics["fixed_b020_iou"]),
            "b020_host_median_iou_mean": float(np.mean(host_metrics["fixed_b020_iou"])),
            "paired_iou_gain_selected_minus_b020_median": median(gains),
            "paired_iou_gain_selected_minus_b020_mean": float(np.mean(gains)),
            "n_hosts_with_positive_iou_gain": int(sum(gain > 0.0 for gain in gains)),
            "selected_completeness_median_across_hosts": median(host_metrics["selected_completeness"]),
            "b020_completeness_median_across_hosts": median(host_metrics["fixed_b020_completeness"]),
            "selected_precision_median_across_hosts": median(host_metrics["selected_precision"]),
            "b020_precision_median_across_hosts": median(host_metrics["fixed_b020_precision"]),
            "n_hosts_with_selected_host_median_iou_below_0p90": int(sum(x < 0.90 for x in selected_host_iou)),
            "n_hosts_with_any_selected_view_iou_below_0p90": int(sum(x < 0.90 for x in selected_min_view_iou)),
            "n_hosts_with_b020_host_median_iou_below_0p90": int(sum(x < 0.90 for x in host_metrics["fixed_b020_iou"])),
            "n_hosts_with_any_b020_view_iou_below_0p90": int(
                sum(min(float(row["fixed_b020_iou"]) for row in by_host[host_id]) < 0.90 for host_id in by_host)
            ),
        }
        summary["bands"][band] = band_summary

    view_csv = HERE / "out_of_fold_view_metrics.csv"
    host_csv = HERE / "out_of_fold_host_metrics.csv"
    folds_csv = HERE / "fold_parameters.csv"
    summary_json = HERE / "summary.json"
    write_csv(view_csv, view_rows)
    write_csv(host_csv, host_rows)
    write_csv(folds_csv, fold_rows)
    summary["analysis_script_sha256"] = sha256(Path(__file__))
    summary["python_version"] = platform.python_version()
    summary["numpy_version"] = np.__version__
    summary_json.write_text(json.dumps(summary, indent=2) + "\n")

    manifest_files = [
        Path(__file__), HERE / "README.md", view_csv, host_csv, folds_csv,
        summary_json, METADATA, CURVES,
    ]
    (HERE / "SHA256SUMS").write_text(
        "".join(f"{sha256(path)}  {os.path.relpath(path, HERE)}\n" for path in manifest_files)
    )
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
