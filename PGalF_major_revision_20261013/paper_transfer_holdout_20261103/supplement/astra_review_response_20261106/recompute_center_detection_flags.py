#!/usr/bin/env python3
"""Recompute central-pixel detection and image-edge flags after distance correction."""
from __future__ import annotations

import argparse
import csv
import json
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
import multiprocessing as mp
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
PAPER = HERE.parents[1]
PROJECT = PAPER.parent
ROOT = PROJECT.parent
SOURCE = PAPER / "scripts" / "source_archive" / "major_revision" / "scripts" / "source"
sys.path.insert(0, str(SOURCE))
import sb_revise as rv  # noqa: E402
import sb_vs_b as sb  # noqa: E402

MANIFEST = PAPER / "supplement" / "sample_manifest.tsv"
PHOTOMETRY = PAPER / "supplement" / "photometry_audit.json"
OUT = HERE / "results" / "detection_flags"


def read_tsv(path: Path):
    with path.open(newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict]):
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    with tmp.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    tmp.replace(path)


def border_touches(region: np.ndarray) -> bool:
    return bool(region[0, :].any() or region[-1, :].any()
                or region[:, 0].any() or region[:, -1].any())


def measure(sample: dict, factors: dict) -> list[dict]:
    step, nstar = sample["step"], int(sample["nstar"])
    h = rv.prepare_halo(step, nstar)
    z = float(h["z"])
    fz = min(factors, key=lambda x: abs(x - z))
    if abs(fz - z) > 2.0e-5:
        raise RuntimeError(f"No audited photometry factor for z={z}")
    flux = np.ascontiguousarray(h["flux"] * factors[fz])
    scale, npix = rv.angular_scale(h["a_snap"], z)
    raw, _pixs = rv.project_raw(h["dx"], h["dy"], h["dz"], flux, scale, npix)
    kernels = [sb.psf_kernel(sb.PIX_AS, sb.FWHM_AS, sb.MOFFAT_BETA, w,
                             sb.KERNEL_AS)[0] for w in sb.WING_A]
    rows = []
    for iview, los in enumerate(rv.LOS):
        for iband, band in enumerate(rv.BAND):
            image = rv.convolve_band(raw[iview][iband], kernels[iband])
            _mu, mask, labels, _ = rv.mu_labels(image, sb.THRESHOLDS[iband])
            cid, fallback, off = rv.component_at(labels, mask, npix // 2, npix // 2)
            region, area = rv.region_of(labels, cid, sb.PIX_AS)
            centre_pixel_detected = bool(labels[npix // 2, npix // 2] > 0)
            empty = not bool(region.any())
            edge = border_touches(region)
            rows.append({
                "step": step, "nstar": nstar, "z": z, "los": los, "band": band,
                "photometry_factor": factors[fz],
                "central_pixel_detected": int(centre_pixel_detected),
                "nearest_bright_pixel_fallback": int(bool(fallback)),
                "empty_contour": int(empty), "touches_image_edge": int(edge),
                "valid_central_isophote": int(centre_pixel_detected and not empty and not edge),
                "area_arcsec2": area, "centre_off_stamp": int(bool(off)),
            })
    return rows


def worker(sample, factors, outdir, index, total):
    step, nstar = sample["step"], int(sample["nstar"])
    print(f"[{index}/{total}] corrected-centre flags {step}_{nstar}", flush=True)
    rows = measure(sample, factors)
    path = outdir / f"{step}_{nstar}.tsv"
    write_tsv(path, rows)
    return rows


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=OUT)
    parser.add_argument("--workers", type=int, default=8)
    args = parser.parse_args()
    outdir = args.output.resolve()
    if outdir.exists() and any(outdir.iterdir()):
        raise FileExistsError(f"Refusing to overwrite {outdir}")
    outdir.mkdir(parents=True, exist_ok=True)
    samples = [r for r in read_tsv(MANIFEST)
               if r["calibration_retained"].lower() == "true"]
    with PHOTOMETRY.open() as fh:
        audit = json.load(fh)
    factors = {round(float(r["z"]), 6): float(r["flux_factor"])
               for r in audit["records"] if "flux_factor" in r}
    all_rows = []
    context = mp.get_context("fork")
    if args.workers <= 1:
        for i, sample in enumerate(samples, 1):
            all_rows.extend(worker(sample, factors, outdir, i, len(samples)))
    else:
        with ProcessPoolExecutor(max_workers=args.workers, mp_context=context) as pool:
            futures = {pool.submit(worker, s, factors, outdir, i, len(samples)): i
                       for i, s in enumerate(samples, 1)}
            finished = {}
            for future in as_completed(futures):
                i = futures[future]
                finished[i] = future.result()
                print(f"[{i}/{len(samples)}] detection flags saved", flush=True)
            for i in sorted(finished):
                all_rows.extend(finished[i])
    write_tsv(outdir / "central_detection_flags.tsv", all_rows)
    print(f"WROTE {len(all_rows)} corrected-distance view flags", flush=True)


if __name__ == "__main__":
    main()
