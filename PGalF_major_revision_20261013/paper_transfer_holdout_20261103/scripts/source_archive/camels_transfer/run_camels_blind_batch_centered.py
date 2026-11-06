#!/usr/bin/env python3
"""Rerun the locked CAMELS sample with every spatial metric centered consistently."""
from __future__ import annotations

import csv
import gc
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed

import measure_camels_halo as worker


TAG = "blind20261030_v4_potential_centered"
MEMBERSHIP_SOURCE = "catalogue_ids"
MAX_WORKERS = 4


def measure_one(task: tuple[int, int]):
    box, group = task
    return worker.measure(
        box, group, save_images=False, tag=TAG,
        membership_source=MEMBERSHIP_SOURCE,
    )


def main() -> None:
    with worker.MANIFEST.open(newline="") as stream:
        manifest = list(csv.DictReader(stream))
    if len(manifest) != 49:
        raise ValueError(f"Expected 49 locked sample rows, found {len(manifest)}")
    by_box: dict[int, list[int]] = defaultdict(list)
    for row in manifest:
        by_box[int(row["cv_box"])].append(int(row["fof_group_index"]))

    for box, groups in by_box.items():
        for group in groups:
            result_dir = worker.RESULTS / f"CV{box:02d}_group{group:05d}_{TAG}"
            if result_dir.exists():
                raise FileExistsError(f"Refusing to overwrite: {result_dir}")

    for box, groups in sorted(by_box.items()):
        snapshot = worker.DATA_ROOT.parent / "raw" / f"CV_{box}" / f"snapshot_{worker.SNAPSHOT:03d}.hdf5"
        catalogue = worker.DATA_ROOT / f"CV_{box}" / f"groups_{worker.SNAPSHOT:03d}.hdf5"
        if not snapshot.is_file() or not catalogue.is_file():
            raise FileNotFoundError(f"Missing CAMELS inputs for CV_{box}: {snapshot} or {catalogue}")
        print(f"Starting CV_{box}: {len(groups)} selected halos on {MAX_WORKERS} workers", flush=True)
        with ProcessPoolExecutor(max_workers=MAX_WORKERS) as pool:
            futures = [pool.submit(measure_one, (box, group)) for group in groups]
            for completed, future in enumerate(as_completed(futures), start=1):
                result = future.result()
                print(f"CV_{box}: completed {completed}/{len(groups)} -> {result}", flush=True)
        gc.collect()
    print(f"Completed {len(manifest)} locked halos with tag {TAG}", flush=True)


if __name__ == "__main__":
    main()
