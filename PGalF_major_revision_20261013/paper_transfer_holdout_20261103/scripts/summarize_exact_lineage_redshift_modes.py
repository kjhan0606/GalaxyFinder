#!/usr/bin/env python3
"""Summarize paired artificial-redshift controls by exact-ID lineage group."""
from __future__ import annotations

import csv
import json
import math
from pathlib import Path

import numpy as np


PAPER = Path(__file__).resolve().parents[1]
SUPP = PAPER / "supplement"
SOURCE = SUPP / "population_redshift_exact_lineage_components_major_20261030.tsv"
SOURCE_METADATA = SUPP / "population_redshift_exact_lineage_metadata_major_20261030.json"
OUTPUT = SUPP / "population_redshift_exact_lineage_paired_modes_major_20261030.tsv"
METADATA = SUPP / "population_redshift_exact_lineage_paired_modes_major_metadata_20261030.json"
DRAWS = 5000
SEED = 20261031
REFERENCE = "observer_full"
COMPARISONS = ("observer_no_tolman_factor", "rest_sed_with_distance")
METRICS = (
    "iou",
    "completeness",
    "precision",
    "non_detection_fraction_los",
    "signed_logL_error",
    "abs_logL_error",
)


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        raise RuntimeError(f"Refusing to write an empty table: {path}")
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter="\t",
                                lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    if OUTPUT.exists() or METADATA.exists():
        raise FileExistsError("Preserve existing results and choose a new output suffix.")
    source_rows = read_tsv(SOURCE)
    rows = [row for row in source_rows if row["b"] == "0.10"]
    lookup: dict[tuple[str, str, str, str], dict[str, str]] = {}
    for row in rows:
        key = (row["z_obs_grid"], row["band"], row["flux_mode"], row["component"])
        if key in lookup:
            raise RuntimeError(f"Duplicate component measurement: {key}")
        lookup[key] = row

    z_values = sorted({row["z_obs_grid"] for row in rows}, key=float)
    bands = sorted({row["band"] for row in rows})
    rng = np.random.default_rng(SEED)
    output_rows: list[dict[str, object]] = []
    for z in z_values:
        for band in bands:
            reference_components = {
                row["component"] for row in rows
                if row["z_obs_grid"] == z and row["band"] == band
                and row["flux_mode"] == REFERENCE
            }
            for comparison in COMPARISONS:
                comparison_components = {
                    row["component"] for row in rows
                    if row["z_obs_grid"] == z and row["band"] == band
                    and row["flux_mode"] == comparison
                }
                common = sorted(reference_components & comparison_components)
                for metric in METRICS:
                    pairs = []
                    for component in common:
                        ref = lookup[(z, band, REFERENCE, component)][metric]
                        alt = lookup[(z, band, comparison, component)][metric]
                        if ref == "" or alt == "":
                            continue
                        ref_value, alt_value = float(ref), float(alt)
                        if math.isfinite(ref_value) and math.isfinite(alt_value):
                            pairs.append((component, alt_value - ref_value))
                    values = np.asarray([value for _, value in pairs], dtype=float)
                    if values.size == 0:
                        continue
                    draw_indices = rng.integers(0, values.size, size=(DRAWS, values.size))
                    means = values[draw_indices].mean(axis=1)
                    lo, hi = np.quantile(means, [0.025, 0.975])
                    output_rows.append({
                        "z_obs": z,
                        "band": band,
                        "b_link": "0.10",
                        "reference_mode": REFERENCE,
                        "comparison_mode": comparison,
                        "metric": metric,
                        "n_exact_id_components": int(values.size),
                        "mean_paired_difference": float(values.mean()),
                        "median_paired_difference": float(np.median(values)),
                        "bootstrap95_low": float(lo),
                        "bootstrap95_high": float(hi),
                    })

    write_tsv(OUTPUT, output_rows)
    source_metadata = json.loads(SOURCE_METADATA.read_text())
    metadata = {
        "description": "Within-component paired contrasts between artificial-redshift image modes at fixed b_link=0.10.",
        "input": SOURCE.name,
        "output": OUTPUT.name,
        "reference_mode": REFERENCE,
        "comparison_modes": list(COMPARISONS),
        "metrics": list(METRICS),
        "resampling_unit": "major-descendant exact-DM-ID connected component; input reduces projections within each output and repeated outputs within each component",
        "components": source_metadata["exact_id_connected_components"],
        "entries": source_metadata["host_output_entries"],
        "exact_id_edges": source_metadata["positive_exact_id_edges_used"],
        "bootstrap_draws": DRAWS,
        "seed": SEED,
        "contrast_sign": "comparison_mode minus observer_full; a negative non_detection_fraction_los difference means fewer views without a detected isophote in the comparison mode",
        "interpretation_limit": source_metadata["interpretation_limit"],
    }
    METADATA.write_text(json.dumps(metadata, indent=2) + "\n")
    print(f"Wrote {len(output_rows)} paired contrasts from {len(z_values)} redshifts, {len(bands)} bands, and {metadata['components']} exact-ID components")
    print(OUTPUT)
    print(METADATA)


if __name__ == "__main__":
    main()
