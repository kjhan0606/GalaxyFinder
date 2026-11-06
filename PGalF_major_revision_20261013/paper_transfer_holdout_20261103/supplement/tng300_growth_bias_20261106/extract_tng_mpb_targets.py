#!/usr/bin/env python3
"""Extract fixed TNG target epochs and audit central status from cached MPBs."""
from __future__ import annotations

import csv
from pathlib import Path

import h5py
import numpy as np

HERE = Path(__file__).resolve().parent
PAPER = HERE.parents[1]
SAMPLE = PAPER / "supplement/tng300_subhalopos_r50_transfer_20261103/tng300_subhalo_pos_40.csv"
TARGET_SNAPSHOTS = (84, 67, 50)
LABELS = {84: "z0p2", 67: "z0p5", 50: "z1p0"}


def read_csv(path: Path):
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream))


def main() -> None:
    if (HERE / "lineage_epochs.csv").exists() or (HERE / "lineage_cohort.csv").exists():
        raise RuntimeError("Refusing to overwrite extracted MPB target tables")
    sample = read_csv(SAMPLE)
    z_by_snap = {int(r["snapshot"]): float(r["redshift"])
                 for r in read_csv(HERE / "snapshot_redshifts.csv")}
    rows = []
    for item in sample:
        rank = int(item["sample_rank"])
        expected_descendant = int(item["central_subhalo_id"])
        path = HERE / f"mpb_rank_{rank:02d}.hdf5"
        if not path.is_file():
            raise FileNotFoundError(path)
        with h5py.File(path, "r") as tree:
            snaps = tree["SnapNum"][:].astype(int)
            for snap in TARGET_SNAPSHOTS:
                matches = np.flatnonzero(snaps == snap)
                if len(matches) != 1:
                    rows.append({
                        "sample_rank": rank,
                        "descendant_group_id": int(item["group_id"]),
                        "descendant_subhalo_id": expected_descendant,
                        "snapshot": snap,
                        "redshift": z_by_snap[snap],
                        "branch_node_present": 0,
                    })
                    continue
                i = int(matches[0])
                subfind_id = int(tree["SubfindID"][i])
                if snap == 84 and subfind_id != expected_descendant:
                    raise RuntimeError(f"Snapshot-84 root mismatch for sample rank {rank}")
                pos = tree["SubhaloPos"][i].astype(float)
                mass_type = tree["SubhaloMassType"][i].astype(float)
                halfmass_type = tree["SubhaloHalfmassRadType"][i].astype(float)
                len_type = tree["SubhaloLenType"][i].astype(int)
                row = {
                    "sample_rank": rank,
                    "descendant_group_id": int(item["group_id"]),
                    "descendant_subhalo_id": expected_descendant,
                    "snapshot": snap,
                    "redshift": z_by_snap[snap],
                    "branch_node_present": 1,
                    "progenitor_subhalo_id": subfind_id,
                    "progenitor_fof_group_id": int(tree["SubhaloGrNr"][i]),
                    "fof_central_subhalo_id": int(tree["GroupFirstSub"][i]),
                    "is_fof_central": int(subfind_id == int(tree["GroupFirstSub"][i])),
                    "subhalo_x_ckpc_h": float(pos[0]),
                    "subhalo_y_ckpc_h": float(pos[1]),
                    "subhalo_z_ckpc_h": float(pos[2]),
                    "stellar_mass_1e10_msun_h": float(mass_type[4]),
                    "stellar_halfmass_radius_ckpc_h": float(halfmass_type[4]),
                    "stellar_particle_count": int(len_type[4]),
                    "subhalo_mass_1e10_msun_h": float(tree["SubhaloMass"][i]),
                    "fof_mass_1e10_msun_h": float(tree["GroupMass"][i]),
                    "fof_mcrit200_1e10_msun_h": float(tree["Group_M_Crit200"][i]),
                    "fof_rcrit200_ckpc_h": float(tree["Group_R_Crit200"][i]),
                }
                rows.append(row)

    with (HERE / "lineage_epochs.csv").open("x", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    by_rank = {}
    for row in rows:
        by_rank.setdefault(row["sample_rank"], []).append(row)
    cohort = []
    for rank, entries in sorted(by_rank.items()):
        present = [r for r in entries if r.get("branch_node_present") == 1]
        cohort.append({
            "sample_rank": rank,
            "descendant_group_id": entries[0]["descendant_group_id"],
            "descendant_subhalo_id": entries[0]["descendant_subhalo_id"],
            "n_target_epochs_present": len(present),
            "n_target_epochs_fof_central": sum(r.get("is_fof_central", 0) == 1 for r in present),
            "all_targets_present": int(len(present) == len(TARGET_SNAPSHOTS)),
            "central_at_all_targets": int(
                len(present) == len(TARGET_SNAPSHOTS)
                and all(r.get("is_fof_central", 0) == 1 for r in present)
            ),
        })
    with (HERE / "lineage_cohort.csv").open("x", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(cohort[0]))
        writer.writeheader()
        writer.writerows(cohort)

    print({
        "n_descendants": len(sample),
        "targets": {LABELS[s]: {"snapshot": s, "z": z_by_snap[s]} for s in TARGET_SNAPSHOTS},
        "branch_present_all_targets": sum(row["all_targets_present"] for row in cohort),
        "central_all_targets": sum(row["central_at_all_targets"] for row in cohort),
        "central_counts": {
            LABELS[s]: sum(r.get("is_fof_central", 0) == 1 for r in rows
                           if r["snapshot"] == s and r.get("branch_node_present") == 1)
            for s in TARGET_SNAPSHOTS
        },
    })


if __name__ == "__main__":
    main()
