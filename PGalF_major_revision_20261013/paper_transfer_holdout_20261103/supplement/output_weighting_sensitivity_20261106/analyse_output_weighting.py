#!/usr/bin/env python3
"""Quantify output-weighting sensitivity of the corrected-distance calibration."""

from __future__ import annotations

import csv
import hashlib
import json
import os
import platform
from collections import defaultdict
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
PAPER = HERE.parents[1]
SUPPLEMENT = PAPER / "supplement"
METADATA = SUPPLEMENT / "halo_metadata.tsv"
CURVES = SUPPLEMENT / "recomputed_curves.tsv"
EXPECTED_ENTRIES = 75
VIEWS = ("x", "y", "z")
BANDS = ("V", "g")


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


def main() -> None:
    metadata = read_tsv(METADATA)
    retained: dict[tuple[str, int], dict[str, str]] = {}
    for row in metadata:
        key = (row["step"], int(row["nstar"]))
        if int(row["interior_V"]) and int(row["interior_g"]):
            retained[key] = row
    if len(retained) != EXPECTED_ENTRIES:
        raise ValueError(f"Expected {EXPECTED_ENTRIES} entries, found {len(retained)}")

    curve_rows = read_tsv(CURVES)
    curves: dict[tuple[str, int, str, str], list[tuple[float, float]]] = defaultdict(list)
    for row in curve_rows:
        key = (row["step"], int(row["nstar"]))
        if row["variant"] == "consistent_distance" and key in retained:
            curves[(*key, row["los"], row["band"])].append(
                (float(row["b"]), float(row["iou"]))
            )

    output_rows: list[dict[str, object]] = []
    summary: dict[str, object] = {
        "analysis": "equal weighting of output-specific median IoU curves",
        "inputs": {
            "halo_metadata.tsv": sha256(METADATA),
            "recomputed_curves.tsv": sha256(CURVES),
        },
        "sample": {
            "n_retained_halo_output_entries": len(retained),
            "n_outputs": len({step for step, _ in retained}),
            "retention": "interior maxima in all three views in both bands",
            "primary_photometry": "consistent_distance, Darwin cosmology",
        },
        "method": {
            "per_entry_curve": "median IoU across x, y, and z at each linking parameter",
            "primary_summary": "median of each entry's median-over-view best parameter",
            "equal_output_curve": "arithmetic mean of the six within-output median-IoU curves",
            "tie_rule": "select the smallest tested parameter at a tied maximum",
            "uncertainty": "descriptive sensitivity; no confidence interval across dependent outputs",
        },
        "bands": {},
    }

    for band in BANDS:
        per_entry_curve: dict[tuple[str, int], np.ndarray] = {}
        per_entry_b: dict[tuple[str, int], float] = {}
        b_grid: np.ndarray | None = None

        for key in retained:
            view_curves = []
            view_best = []
            for los in VIEWS:
                rows = sorted(curves[(*key, los, band)])
                current_b = np.asarray([item[0] for item in rows], dtype=float)
                current_iou = np.asarray([item[1] for item in rows], dtype=float)
                if b_grid is None:
                    b_grid = current_b
                if not np.array_equal(current_b, b_grid):
                    raise ValueError(f"Inconsistent parameter grid for {key}, {los}, {band}")
                view_curves.append(current_iou)
                view_best.append(float(current_b[int(np.argmax(current_iou))]))
            if len(view_curves) != len(VIEWS):
                raise ValueError(f"Missing projections for {key}, {band}")
            per_entry_curve[key] = np.median(np.stack(view_curves), axis=0)
            per_entry_b[key] = float(np.median(view_best))

        if b_grid is None:
            raise ValueError(f"No corrected curves for band {band}")

        steps = sorted({key[0] for key in retained})
        output_median_curves = []
        output_mean_curves = []
        output_level_best: dict[str, float] = {}
        output_level_rows = []
        for step in steps:
            keys = [key for key in retained if key[0] == step]
            matrix = np.stack([per_entry_curve[key] for key in keys])
            median_curve = np.median(matrix, axis=0)
            mean_curve = np.mean(matrix, axis=0)
            output_median_curves.append(median_curve)
            output_mean_curves.append(mean_curve)
            median_best = float(b_grid[int(np.argmax(median_curve))])
            mean_best = float(b_grid[int(np.argmax(mean_curve))])
            output_level_best[step] = median_best
            z = float(retained[keys[0]]["z"])
            output_level_rows.append(
                {
                    "band": band,
                    "step": step,
                    "redshift": z,
                    "n_entries": len(keys),
                    "median_entry_best_b": float(np.median([per_entry_b[key] for key in keys])),
                    "within_output_median_curve_best_b": median_best,
                    "within_output_mean_curve_best_b": mean_best,
                    "within_output_median_curve_max_iou": float(np.max(median_curve)),
                }
            )

        equal_output_median_curve = np.mean(np.stack(output_median_curves), axis=0)
        equal_output_mean_curve = np.mean(np.stack(output_mean_curves), axis=0)
        primary_b = float(np.median(list(per_entry_b.values())))
        equal_output_b = float(b_grid[int(np.argmax(equal_output_median_curve))])
        equal_output_mean_b = float(b_grid[int(np.argmax(equal_output_mean_curve))])

        output_rows.extend(output_level_rows)
        summary["bands"][band] = {
            "n_entries": len(retained),
            "tested_b_grid": b_grid.tolist(),
            "entry_weighted_median_of_entry_best_b": primary_b,
            "equal_output_weighted_best_b": equal_output_b,
            "equal_output_weighted_best_b_using_output_mean_curves": equal_output_mean_b,
            "median_of_six_output_specific_median_curve_best_b": float(
                np.median(list(output_level_best.values()))
            ),
            "equal_output_curve_max_iou": float(np.max(equal_output_median_curve)),
            "equal_output_curve_iou_at_entry_weighted_b": float(
                equal_output_median_curve[int(np.flatnonzero(b_grid == primary_b)[0])]
            ),
            "equal_output_curve_iou_at_equal_output_b": float(np.max(equal_output_median_curve)),
        }

    output_csv = HERE / "output_weighting_by_output.csv"
    summary_json = HERE / "output_weighting_summary.json"
    write_csv(output_csv, output_rows)
    summary["analysis_script_sha256"] = sha256(Path(__file__))
    summary["python_version"] = platform.python_version()
    summary["numpy_version"] = np.__version__
    summary_json.write_text(json.dumps(summary, indent=2) + "\n")

    manifest_files = [Path(__file__), HERE / "README.md", output_csv, summary_json, METADATA, CURVES]
    (HERE / "SHA256SUMS").write_text(
        "".join(f"{sha256(path)}  {os.path.relpath(path, HERE)}\n" for path in manifest_files)
    )
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
