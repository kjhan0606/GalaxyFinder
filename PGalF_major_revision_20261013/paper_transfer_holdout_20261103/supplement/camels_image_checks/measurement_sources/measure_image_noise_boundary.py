#!/usr/bin/env python3
"""Stronger fixed-calibration noise stress with rms equal to contour pixel flux."""
from copy import deepcopy
import csv
import gc
import hashlib
import json
from concurrent.futures import ProcessPoolExecutor,as_completed
from pathlib import Path
import measure_image_sensitivity as worker

ROOT=Path(__file__).resolve().parent
OUT=ROOT/"results/image_noise_boundary_20261030"
BASE_CONFIG=deepcopy(worker.CONFIG)


def scenarios():
    return [("baseline",.2,0.,False,None,None,0)]+[
        ("noise_1sigma",.2,0.,False,None,1.,seed) for seed in range(3)]


def configure():
    config=deepcopy(BASE_CONFIG)
    config.update(active_conditions=["baseline","noise_1sigma"],threshold_offsets_mag=[0.],
                  coarse_pixel_arcsec=None,no_wing="not varied",sky_residual_mag_arcsec2=None,
                  noise_threshold_to_pixel_rms=[1.],
                  test_scope="Explicit stronger-noise extension; no linking-parameter retuning")
    worker.CONFIG=config
    worker.CONFIG_SHA=hashlib.sha256(json.dumps(config,sort_keys=True).encode()).hexdigest()
    worker.OUT=OUT
    worker.scenarios=scenarios


def measure_one(task):
    configure()
    return worker.measure(task)


def main():
    configure()
    OUT.mkdir(parents=True,exist_ok=True)
    config=OUT/"configuration.json"
    if config.exists() and json.loads(config.read_text())["config_sha256"]!=worker.CONFIG_SHA:
        raise ValueError("Frozen strong-noise configuration differs")
    config.write_text(json.dumps({"config_sha256":worker.CONFIG_SHA,**worker.CONFIG},indent=2)+"\n")
    manifest=list(csv.DictReader(worker.base.MANIFEST.open()))
    if len(manifest)!=49:raise ValueError("Locked manifest changed")
    tasks=[(int(r["cv_box"]),int(r["fof_group_index"])) for r in manifest]
    print("pilot",measure_one(tasks[0]),flush=True)
    for box in sorted(set(b for b,_ in tasks)):
        with ProcessPoolExecutor(max_workers=2) as pool:
            futures=[pool.submit(measure_one,t) for t in tasks if t[0]==box]
            for future in as_completed(futures):print("completed",future.result(),flush=True)
        worker.base.clear_snapshot_id_cache()
        gc.collect()


if __name__=="__main__":main()
