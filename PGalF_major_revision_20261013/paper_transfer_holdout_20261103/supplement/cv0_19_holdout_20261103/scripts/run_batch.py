#!/usr/bin/env python3
"""Measure all frozen CV0--19 hosts using the score-unseen holdout setup."""
from __future__ import annotations

import csv
import gc
import json
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import measure_camels_halo_holdout as worker


ROOT = Path(__file__).resolve().parent
TAG = "cv019_holdout_20261103"
MEMBERSHIP_SOURCE = "catalogue_ids"
MAX_WORKERS = 4
worker.MANIFEST = ROOT / "sample_manifest.csv"
worker.RESULTS = ROOT / "results"


def measure_one(task: tuple[int, int]):
    box, group = task
    return worker.measure(
        box, group, save_images=False, tag=TAG,
        membership_source=MEMBERSHIP_SOURCE,
    )


def main() -> None:
    with worker.MANIFEST.open(newline="") as stream:
        manifest = list(csv.DictReader(stream))
    if len(manifest) != 130:
        raise ValueError(f"Expected 130 frozen hosts, found {len(manifest)}")
    by_box: dict[int, list[int]] = defaultdict(list)
    for row in manifest:
        by_box[int(row["cv_box"])].append(int(row["fof_group_index"]))
    if sorted(by_box) != list(range(20)):
        raise ValueError(f"Expected CV0--CV19; found boxes {sorted(by_box)}")

    for box, groups in by_box.items():
        snapshot = Path("/gpfs/kjhan/CAMELS/Swift-EAGLE/L25n256/raw") / f"CV_{box}" / f"snapshot_{worker.SNAPSHOT:03d}.hdf5"
        catalog = worker.DATA_ROOT / f"CV_{box}" / f"groups_{worker.SNAPSHOT:03d}.hdf5"
        if not snapshot.is_file() or not catalog.is_file():
            raise FileNotFoundError(f"Missing CV_{box} inputs: {snapshot} or {catalog}")
        for group in groups:
            result = worker.RESULTS / f"CV{box:02d}_group{group:05d}_{TAG}"
            if result.exists():
                metrics = result / "metrics.csv"
                metadata = result / "metadata.json"
                if not metrics.is_file() or not metadata.is_file():
                    raise FileExistsError(f"Refusing to overwrite incomplete result: {result}")
                with metrics.open(newline="") as stream:
                    reader = csv.DictReader(stream)
                    count = sum(1 for _ in reader)
                    if "target_touches_image_edge" not in (reader.fieldnames or []):
                        raise ValueError(f"Incompatible existing metrics file: {metrics}")
                with metadata.open() as stream:
                    saved = json.load(stream)
                if count != 288 or int(saved["cv_box"]) != box or int(saved["fof_group_index"]) != group:
                    raise ValueError(f"Existing result is incomplete or mismatched: {result}")

    print(f"Frozen manifest: {len(manifest)} hosts in {len(by_box)} boxes", flush=True)
    print(f"Membership source: {MEMBERSHIP_SOURCE}; workers per box: {MAX_WORKERS}", flush=True)
    completed_total = 0
    for box, groups in sorted(by_box.items()):
        tasks = []
        for group in groups:
            result = worker.RESULTS / f"CV{box:02d}_group{group:05d}_{TAG}"
            if result.exists():
                completed_total += 1
            else:
                tasks.append((box, group))
        if not tasks:
            print(f"Skipping CV_{box}: all {len(groups)} hosts already complete", flush=True)
            continue
        print(f"Starting CV_{box}: {len(tasks)} hosts ({len(groups) - len(tasks)} already complete)", flush=True)
        with ProcessPoolExecutor(max_workers=MAX_WORKERS) as pool:
            futures = [pool.submit(measure_one, task) for task in tasks]
            for completed, future in enumerate(as_completed(futures), start=1):
                print(f"CV_{box}: completed {completed}/{len(tasks)} -> {future.result()}", flush=True)
                completed_total += 1
        gc.collect()
    print(f"Complete: {completed_total}/{len(manifest)} frozen hosts", flush=True)


if __name__ == "__main__":
    main()
