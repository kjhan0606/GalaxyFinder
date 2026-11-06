#!/usr/bin/env python3
"""Leave-one-output check using the cosmology-corrected calibration curves."""
from __future__ import annotations

import csv
import json
from collections import defaultdict
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
SUPP = ROOT / "supplement"
B_GRID = np.array([0.02, 0.03, 0.05, 0.08, 0.10, 0.12, 0.15, 0.18, 0.20])
BOX = 44.171267


def read_tsv(path: Path):
    with path.open() as stream:
        return list(csv.DictReader(
            (line for line in stream if line.strip() and not line.startswith("#")),
            delimiter="\t",
        ))


def main() -> None:
    metadata = read_tsv(SUPP / "halo_metadata.tsv")
    positions = {
        (row["step"], int(row["nstar"])): np.array(
            [float(row[key]) for key in ("x_hinv_cMpc", "y_hinv_cMpc", "z_hinv_cMpc")]
        )
        for row in read_tsv(SUPP / "halo_positions.tsv")
    }
    keys = [(row["step"], int(row["nstar"])) for row in metadata]
    retained = {
        (row["step"], int(row["nstar"])):
        bool(int(row["interior_V"]) and int(row["interior_g"]))
        for row in metadata
    }
    curves = defaultdict(list)
    for row in read_tsv(SUPP / "recomputed_curves.tsv"):
        if row["variant"] != "consistent_distance":
            continue
        key = (row["step"], int(row["nstar"]), row["los"], row["band"])
        curves[key].append(row)

    halo_curves = {}
    for key in keys:
        for band in ("V", "g"):
            views = []
            for los in "xyz":
                rows = sorted(curves[(*key, los, band)], key=lambda row: float(row["b"]))
                if len(rows) != len(B_GRID):
                    raise ValueError(f"Incomplete corrected curve for {key} {los} {band}")
                if not np.allclose([float(row["b"]) for row in rows], B_GRID):
                    raise ValueError(f"Unexpected linking grid for {key} {los} {band}")
                views.append(np.array([float(row["iou"]) for row in rows]))
            halo_curves[(*key, band)] = np.median(views, axis=0)

    folds = []
    predictions = []
    for held_step in sorted({key[0] for key in keys}):
        held = [key for key in keys if key[0] == held_step]
        training = []
        nearby_removed = 0
        for key in keys:
            if key[0] == held_step or not retained[key]:
                continue
            displacement = np.array([positions[key] - positions[h] for h in held])
            displacement -= BOX * np.rint(displacement / BOX)
            if np.min(np.linalg.norm(displacement, axis=1)) < 0.5:
                nearby_removed += 1
            else:
                training.append(key)
        for band in ("V", "g"):
            objective = np.median([halo_curves[(*key, band)] for key in training], axis=0)
            selected_index = int(np.argmax(objective))
            selected_b = float(B_GRID[selected_index])
            folds.append({
                "held_out_step": held_step, "band": band,
                "n_train": len(training), "n_neighbours_removed": nearby_removed,
                "n_held": len(held), "selected_b": selected_b,
            })
            for key in held:
                curve = halo_curves[(*key, band)]
                predictions.append({
                    "step": key[0], "nstar": key[1], "band": band,
                    "interior_in_both_bands": int(retained[key]),
                    "selected_b": selected_b,
                    "held_out_iou": float(curve[selected_index]),
                    "fixed_b010_iou": float(curve[int(np.flatnonzero(B_GRID == 0.10)[0])]),
                    "fixed_b020_iou": float(curve[-1]),
                })

    summary = {}
    for band in ("V", "g"):
        rows = [row for row in predictions
                if row["band"] == band and row["interior_in_both_bands"]]
        summary[band] = {
            "n_entries": len(rows),
            "held_out_iou_median": float(np.median([row["held_out_iou"] for row in rows])),
            "gain_over_b020_median": float(np.median([
                row["held_out_iou"] - row["fixed_b020_iou"] for row in rows
            ])),
            "selected_b_values": sorted({row["selected_b"] for row in folds if row["band"] == band}),
        }

    for filename, rows in (("validation_corrected_folds.tsv", folds),
                           ("validation_corrected_predictions.tsv", predictions)):
        with (SUPP / filename).open("w") as stream:
            writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter="\t",
                                    lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)
    (SUPP / "validation_corrected_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n"
    )
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
