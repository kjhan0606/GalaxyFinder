#!/usr/bin/env python3
"""Audit CAMELS catalogue-ID versus snapshot FoF membership diagnostics."""
from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parent
TAG = "cv019_holdout_20261103"


def main() -> None:
    with (ROOT / "sample_manifest.csv").open(newline="") as stream:
        hosts = list(csv.DictReader(stream))
    rows = []
    for host in hosts:
        box = int(host["cv_box"])
        group = int(host["fof_group_index"])
        path = ROOT / "results" / f"CV{box:02d}_group{group:05d}_{TAG}" / "metadata.json"
        meta = json.loads(path.read_text())
        cat_count = int(meta["catalogue_fof_star_count"])
        raw_count = int(meta["raw_fof_star_count"])
        intersection = int(meta["catalogue_group_star_ids_matching_raw_fof"])
        union = cat_count + raw_count - intersection
        satellite_count = int(meta["resolved_satellite_catalogue_star_count"])
        satellite_intersection = int(meta["resolved_satellite_ids_matching_raw_fof"])
        rows.append({
            "cv_box": box,
            "fof_group_index": group,
            "fof_mass_1e10_msun_per_h": float(host["fof_mass_1e10_msun_per_h"]),
            "catalogue_fof_star_count": cat_count,
            "raw_fof_star_count": raw_count,
            "group_id_intersection": intersection,
            "group_id_jaccard": intersection / union if union else float("nan"),
            "catalogue_group_ids_matching_snapshot": int(meta["catalogue_group_star_ids_matching_snapshot"]),
            "catalogue_group_id_resolution_exact": int(meta["catalogue_group_star_ids_matching_snapshot"]) == cat_count,
            "central_star_count": int(meta["central_catalogue_id_count"]),
            "central_ids_matching_raw_fof": int(meta["central_catalogue_ids_matching_raw_fof"]),
            "central_id_match_fraction": float(meta["central_catalogue_id_match_fraction"]),
            "central_ids_unassigned_to_fof": int(meta["central_catalogue_ids_unassigned_to_fof"]),
            "central_ids_assigned_to_other_fof": int(meta["central_catalogue_ids_assigned_to_other_fof"]),
            "resolved_satellite_subhalo_count": int(meta["resolved_satellite_subhalo_count"]),
            "resolved_satellite_catalogue_star_count": satellite_count,
            "resolved_satellite_ids_matching_raw_fof": satellite_intersection,
            "resolved_satellite_id_match_fraction": satellite_intersection / satellite_count if satellite_count else float("nan"),
            "input_star_count_after_satellite_mask": int(meta["input_star_count_after_satellite_mask"]),
        })
    output = ROOT / "membership_audit.csv"
    with output.open("x", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    jaccard = np.asarray([row["group_id_jaccard"] for row in rows], dtype=float)
    central = np.asarray([row["central_id_match_fraction"] for row in rows], dtype=float)
    satellite = np.asarray([row["resolved_satellite_id_match_fraction"] for row in rows], dtype=float)
    report = {
        "n_hosts": len(rows),
        "catalogue_group_id_resolution_failures": sum(not row["catalogue_group_id_resolution_exact"] for row in rows),
        "group_membership_jaccard_median": float(np.median(jaccard)),
        "group_membership_jaccard_p16": float(np.quantile(jaccard, 0.16)),
        "group_membership_jaccard_p84": float(np.quantile(jaccard, 0.84)),
        "group_membership_jaccard_min": float(np.min(jaccard)),
        "group_membership_jaccard_below_0p98": int(np.sum(jaccard < 0.98)),
        "central_id_match_fraction_median": float(np.median(central)),
        "central_id_match_fraction_min": float(np.min(central)),
        "central_id_match_fraction_below_0p99": int(np.sum(central < 0.99)),
        "resolved_satellite_id_match_fraction_median": float(np.nanmedian(satellite)),
        "hosts_with_zero_resolved_satellite_ids_in_raw_fof": [
            {"cv_box": row["cv_box"], "fof_group_index": row["fof_group_index"]}
            for row in rows if row["resolved_satellite_catalogue_star_count"] and not row["resolved_satellite_ids_matching_raw_fof"]
        ],
        "metadata_source": "per-host metadata.json from the locked CV0--CV19 run",
    }
    (ROOT / "membership_audit_summary.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
