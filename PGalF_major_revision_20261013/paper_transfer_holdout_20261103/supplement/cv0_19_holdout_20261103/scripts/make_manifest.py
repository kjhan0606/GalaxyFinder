#!/usr/bin/env python3
"""Create the frozen CV0--19 Swift-EAGLE holdout host manifest."""
from __future__ import annotations

import csv
from pathlib import Path

import h5py
import numpy as np


ROOT = Path(__file__).resolve().parent
DATA_ROOT = Path("/gpfs/kjhan/CAMELS/Swift-EAGLE/L25n256/CV")
OUTPUT = ROOT / "sample_manifest.csv"
SNAPSHOT = 90
MASS_MIN = 1.0e3
MASS_MAX = 1.0e4
MIN_CENTRAL_STARS = 100
FIELDS = [
    "cv_box", "snapshot", "redshift", "fof_group_index", "fof_group_id_raw",
    "fof_mass_1e10_msun_per_h", "fof_star_count", "subhalo_count",
    "central_subhalo_index", "central_star_count",
    "central_stellar_mass_1e10_msun_per_h",
]


def main() -> None:
    if OUTPUT.exists():
        raise FileExistsError(f"Refusing to overwrite {OUTPUT}")
    rows: list[dict[str, int | float]] = []
    counts: dict[int, int] = {}
    for cv_box in range(20):
        source = DATA_ROOT / f"CV_{cv_box}" / f"groups_{SNAPSHOT:03d}.hdf5"
        if not source.is_file():
            raise FileNotFoundError(source)
        with h5py.File(source, "r") as catalog:
            redshift = float(catalog["Header"].attrs["Redshift"])
            group_mass = catalog["Group/GroupMass"][:]
            group_star_count = catalog["Group/GroupLenType"][:, 4]
            first_subhalo = catalog["Group/GroupFirstSub"][:]
            n_subhalos = catalog["Group/GroupNsubs"][:]
            subhalo_len_type = catalog["Subhalo/SubhaloLenType"][:]
            subhalo_mass_type = catalog["Subhalo/SubhaloMassType"][:]
        selected_groups = np.flatnonzero((group_mass >= MASS_MIN) & (group_mass < MASS_MAX))
        box_rows = 0
        for group_index in selected_groups:
            first = int(first_subhalo[group_index])
            count = int(n_subhalos[group_index])
            if count < 1:
                continue
            subhalos = np.arange(first, first + count, dtype=np.int64)
            star_counts = subhalo_len_type[subhalos, 4]
            central_offset = int(np.argmax(star_counts))
            central_index = int(subhalos[central_offset])
            central_stars = int(star_counts[central_offset])
            if central_stars < MIN_CENTRAL_STARS:
                continue
            rows.append({
                "cv_box": cv_box,
                "snapshot": SNAPSHOT,
                "redshift": redshift,
                "fof_group_index": int(group_index),
                "fof_group_id_raw": int(group_index) + 1,
                "fof_mass_1e10_msun_per_h": float(group_mass[group_index]),
                "fof_star_count": int(group_star_count[group_index]),
                "subhalo_count": count,
                "central_subhalo_index": central_index,
                "central_star_count": central_stars,
                "central_stellar_mass_1e10_msun_per_h": float(subhalo_mass_type[central_index, 4]),
            })
            box_rows += 1
        counts[cv_box] = box_rows
    rows.sort(key=lambda row: (int(row["cv_box"]), int(row["fof_group_index"])))
    if len(rows) != 130 or len(counts) != 20:
        raise RuntimeError(f"Frozen selection changed: {len(rows)} hosts in {len(counts)} boxes")
    with OUTPUT.open("x", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=FIELDS)
        writer.writeheader()
        writer.writerows(rows)
    print(f"Wrote {len(rows)} hosts in {len(counts)} CV boxes to {OUTPUT}")
    print("Per-box counts:", ", ".join(f"CV_{box}={count}" for box, count in counts.items()))


if __name__ == "__main__":
    main()
