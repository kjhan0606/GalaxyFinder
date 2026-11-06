#!/usr/bin/env python3
"""Summarize catalogued stellar growth for the locked common-central trees."""
from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
H = 0.6774
SNAPS = (84, 67, 50)
SEED = 8172026
N_BOOT = 20000


def read(path: Path):
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream))


def bootstrap_summary(values: np.ndarray, rng: np.random.Generator):
    values = np.asarray(values, dtype=float)
    samples = rng.integers(0, values.size, size=(N_BOOT, values.size))
    draws = values[samples]
    output = {}
    for name, statistic in (("mean", np.mean), ("median", np.median)):
        point = float(statistic(values))
        interval = np.quantile(statistic(draws, axis=1), [0.025, 0.975])
        output[name] = {"estimate": point, "ci95": [float(interval[0]), float(interval[1])],
                        "unit": "dex"}
    return output


def main():
    cohort = {int(row["sample_rank"]) for row in read(HERE / "lineage_cohort.csv")
              if int(row["central_at_all_targets"]) == 1}
    rows = [row for row in read(HERE / "lineage_epochs.csv")
            if int(row["sample_rank"]) in cohort and int(row["snapshot"]) in SNAPS]
    by_rank = {rank: {} for rank in sorted(cohort)}
    for row in rows:
        rank = int(row["sample_rank"])
        snap = int(row["snapshot"])
        if int(row["is_fof_central"]) != 1:
            raise RuntimeError(f"Non-central branch node in locked sample: rank={rank}, snap={snap}")
        by_rank[rank][snap] = row
    if any(set(values) != set(SNAPS) for values in by_rank.values()):
        raise RuntimeError("A common-central tree is missing a target epoch")

    metrics = []
    for rank, values in by_rank.items():
        row84, row67, row50 = (values[snap] for snap in SNAPS)
        m84 = float(row84["stellar_mass_1e10_msun_h"])
        m67 = float(row67["stellar_mass_1e10_msun_h"])
        m50 = float(row50["stellar_mass_1e10_msun_h"])
        r84 = float(row84["stellar_halfmass_radius_ckpc_h"]) / (H * (1.0 + float(row84["redshift"])))
        r67 = float(row67["stellar_halfmass_radius_ckpc_h"]) / (H * (1.0 + float(row67["redshift"])))
        r50 = float(row50["stellar_halfmass_radius_ckpc_h"]) / (H * (1.0 + float(row50["redshift"])))
        metrics.append({
            "sample_rank": rank,
            "descendant_group_id": int(row84["descendant_group_id"]),
            "descendant_subhalo_id": int(row84["descendant_subhalo_id"]),
            "snap84_redshift": float(row84["redshift"]),
            "snap67_redshift": float(row67["redshift"]),
            "snap50_redshift": float(row50["redshift"]),
            "snap84_stellar_mass_msun": m84 * 1.0e10 / H,
            "snap67_stellar_mass_msun": m67 * 1.0e10 / H,
            "snap50_stellar_mass_msun": m50 * 1.0e10 / H,
            "snap84_stellar_halfmass_radius_pkpc": r84,
            "snap67_stellar_halfmass_radius_pkpc": r67,
            "snap50_stellar_halfmass_radius_pkpc": r50,
            "delta_log10_stellar_mass_z1_to_z0p2_dex": float(np.log10(m84 / m50)),
            "delta_log10_stellar_mass_z0p5_to_z0p2_dex": float(np.log10(m84 / m67)),
            "delta_log10_stellar_mass_z1_to_z0p5_dex": float(np.log10(m67 / m50)),
            "delta_log10_stellar_halfmass_radius_z1_to_z0p2_dex": float(np.log10(r84 / r50)),
            "snap84_stellar_particle_count": int(row84["stellar_particle_count"]),
            "snap67_stellar_particle_count": int(row67["stellar_particle_count"]),
            "snap50_stellar_particle_count": int(row50["stellar_particle_count"]),
        })

    output_csv = HERE / "tree_catalog_growth_metrics.csv"
    output_json = HERE / "tree_catalog_growth_summary.json"
    if output_csv.exists() or output_json.exists():
        raise FileExistsError("Refusing to overwrite existing catalog-growth summaries")
    with output_csv.open("x", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(metrics[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(metrics)

    rng = np.random.default_rng(SEED)
    summaries = {}
    for key in ("delta_log10_stellar_mass_z1_to_z0p2_dex",
                "delta_log10_stellar_mass_z0p5_to_z0p2_dex",
                "delta_log10_stellar_mass_z1_to_z0p5_dex",
                "delta_log10_stellar_halfmass_radius_z1_to_z0p2_dex"):
        summaries[key] = bootstrap_summary(np.array([m[key] for m in metrics]), rng)
    result = {
        "description": "SubLink main-progenitor catalog growth; common FoF-central cohort at all three snapshots",
        "simulation": "TNG300-1",
        "snapshots": {str(s): {"redshift": float(by_rank[next(iter(by_rank))][s]["redshift"]),
                                "n_trees": len(metrics)} for s in SNAPS},
        "n_independent_descendant_trees": len(metrics),
        "unit_of_resampling": "descendant tree (one resampled value per final host)",
        "bootstrap_seed": SEED,
        "bootstrap_replicates": N_BOOT,
        "stellar_mass_definition": "SubhaloMassType[4] * 1e10 / h Msun; catalogued bound stellar mass",
        "size_definition": "3D SubhaloHalfmassRadType[4] converted from ckpc/h to pkpc",
        "growth_summaries": summaries,
        "per_tree_csv": output_csv.name,
    }
    output_json.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({"n_trees": len(metrics), "growth_summaries": summaries}, indent=2))


if __name__ == "__main__":
    main()
