#!/usr/bin/env python3
"""Recompute corrected-distance PSF scores and boundary-width sensitivity.

The corrected luminosity-distance factor is applied before the image is
thresholded. Existing coarse-grid FoF memberships are reused; the six
intermediate linking lengths are calculated from the Darwin particle data.
All outputs are written below this analysis directory.
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
import multiprocessing as mp
from pathlib import Path

import numpy as np
from scipy import ndimage
from scipy.signal import fftconvolve

HERE = Path(__file__).resolve().parent
PAPER = HERE.parents[1]
PROJECT = PAPER.parent
ROOT = PROJECT.parent
SOURCE = PAPER / "scripts" / "source_archive" / "major_revision" / "scripts" / "source"
sys.path.insert(0, str(SOURCE))
import sb_revise as rv  # noqa: E402
import sb_vs_b as sb  # noqa: E402

RESULTS = HERE / "results"
MEMBERSHIP_ROOT = PROJECT / "results"
MANIFEST = PAPER / "supplement" / "sample_manifest.tsv"
CURVES = PAPER / "supplement" / "recomputed_curves.tsv"
BOUNDARY = PAPER / "supplement" / "boundary_ring_corrected.tsv"
PHOTOMETRY = PAPER / "supplement" / "photometry_audit.json"
B_FINE = np.asarray((0.02, 0.03, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10,
                     0.11, 0.12, 0.13, 0.14, 0.15, 0.18, 0.20), dtype=float)
BOUNDARY_DISTANCE_THRESHOLDS_AS = (0.2, 0.3, 0.5, 1.0, 2.0)
BOUNDARY_PIXEL = 0.2


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict]) -> None:
    if not rows:
        raise RuntimeError(f"No rows to write to {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    with tmp.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    tmp.replace(path)


def score_fractional(w: np.ndarray, keep: np.ndarray, q: np.ndarray) -> tuple[float, float, float]:
    """IoU, completeness, precision with fractional image-aperture weights."""
    good = np.isfinite(w) & (w > 0.0)
    wk = np.where(good & keep, w, 0.0)
    wq = np.where(good, w * q, 0.0)
    wbq = np.where(good & keep, w * q, 0.0)
    lg = float(wk.sum())
    li = float(wq.sum())
    lb = float(wbq.sum())
    den = lg + li - lb
    return (lb / den if den > 0.0 else np.nan,
            lb / li if li > 0.0 else np.nan,
            lb / lg if lg > 0.0 else np.nan)


def score_binary(w: np.ndarray, keep: np.ndarray, inside: np.ndarray) -> tuple[float, float, float]:
    good = np.isfinite(w) & (w > 0.0)
    lg = float(w[good & keep].sum())
    li = float(w[good & inside].sum())
    lb = float(w[good & keep & inside].sum())
    den = lg + li - lb
    return (lb / den if den > 0.0 else np.nan,
            lb / li if li > 0.0 else np.nan,
            lb / lg if lg > 0.0 else np.nan)


def score_boundary(w: np.ndarray, keep: np.ndarray, inside: np.ndarray,
                   band: np.ndarray) -> tuple[float, float, float]:
    return score_binary(w[band], keep[band], inside[band])


def signed_distance_arcsec(region: np.ndarray) -> np.ndarray:
    """Pixel-centre distance: negative inside and positive outside."""
    inside = ndimage.distance_transform_edt(region)
    outside = ndimage.distance_transform_edt(~region)
    return (outside - inside) * BOUNDARY_PIXEL


def load_references():
    corrected = {}
    for row in read_tsv(CURVES):
        if row["variant"] != "consistent_distance":
            continue
        key = (row["step"], int(row["nstar"]), row["los"], row["band"],
               round(float(row["b"]), 2))
        corrected[key] = {name: float(row[name]) for name in
                          ("iou", "completeness", "precision")}
    primary = {}
    for row in read_tsv(BOUNDARY):
        key = (row["step"], int(row["nstar"]), row["los"], row["band"])
        primary[key] = float(row["b_best"])
    with PHOTOMETRY.open() as fh:
        audit = json.load(fh)
    factors = {round(float(r["z"]), 6): float(r["flux_factor"])
               for r in audit["records"] if "flux_factor" in r}
    return corrected, primary, factors


def fine_memberships(step: str, nstar: int, mass: np.ndarray,
                     coords: tuple[np.ndarray, np.ndarray, np.ndarray]) -> dict[float, np.ndarray]:
    cache_path = MEMBERSHIP_ROOT / f"{step}_{nstar}" / "membership.npz"
    if not cache_path.is_file():
        raise FileNotFoundError(cache_path)
    with np.load(cache_path) as cache:
        if int(cache["n_stars"]) != len(mass):
            raise RuntimeError(f"Membership cache length mismatch for {step}_{nstar}")
        seed = int(cache["seed"])
        keeps = {}
        for b in B_FINE:
            key = f"b{b:.2f}"
            if key in cache.files:
                keep = np.unpackbits(cache[key], count=len(mass)).astype(bool)
                keeps[float(b)] = keep
        if len(keeps) != 9:
            raise RuntimeError(f"Expected 9 cached coarse-grid groups for {step}_{nstar}")
    # The locked cache stores the density seed used to create the published
    # memberships. Reuse it directly; recomputing the 32-neighbour density
    # for multi-million-particle catalogues would duplicate a costly step.
    if not 0 <= seed < len(mass):
        raise RuntimeError(f"Cached density seed is outside {step}_{nstar}")
    for b in B_FINE:
        b = float(b)
        if b in keeps:
            continue
        root_labels = sb.fof_largest(*coords, sb.link02(mass, b))
        keeps[b] = sb.peak_mask(root_labels, *coords, seed)
        del root_labels
        print(f"    computed inserted b={b:.2f}; group stars={int(keeps[b].sum())}", flush=True)
    return keeps


def measure_halo(step: str, nstar: int, corrected_ref: dict,
                 primary_b: dict, factors: dict) -> list[dict]:
    h = rv.prepare_halo(step, nstar)
    z = float(h["z"])
    factor_key = min(factors, key=lambda value: abs(value - z))
    if abs(factor_key - z) > 2.0e-5:
        raise RuntimeError(f"No audited photometry factor for z={z}")
    factor = factors[factor_key]
    flux = np.ascontiguousarray(h["flux"] * factor)
    coords = (h["dx"], h["dy"], h["dz"])
    keeps = fine_memberships(step, nstar, h["mass"], coords)
    scale, npix = rv.angular_scale(h["a_snap"], z)
    raw, pixs = rv.project_raw(*coords, flux, scale, npix)
    kernels = [sb.psf_kernel(sb.PIX_AS, sb.FWHM_AS, sb.MOFFAT_BETA,
                             wing, sb.KERNEL_AS)[0] for wing in sb.WING_A]
    rows: list[dict] = []

    for iview, los in enumerate(rv.LOS):
        ix, iy, on = pixs[iview]
        for iband, band in enumerate(rv.BAND):
            image = rv.convolve_band(raw[iview][iband], kernels[iband])
            _mu, mask, labels, _ = rv.mu_labels(image, sb.THRESHOLDS[iband])
            cid, _fallback, _off = rv.component_at(labels, mask, npix // 2, npix // 2)
            region, area = rv.region_of(labels, cid, sb.PIX_AS)
            inside = rv.sample_inside(region, ix, iy, on)

            if region.any():
                q_map = fftconvolve(region.astype(np.float64), kernels[iband], mode="same")
                q = np.zeros(len(flux), dtype=np.float64)
                q[on] = q_map[iy[on], ix[on]]
                q = np.clip(q, 0.0, 1.0)
                aperture_from_particles = float(np.sum(flux[:, iband] * q))
                aperture_from_image = float(image[region].sum())
                scale_check = max(abs(aperture_from_image), abs(aperture_from_particles), 1.0e-300)
                if abs(aperture_from_particles - aperture_from_image) > 2.0e-4 * scale_check:
                    raise RuntimeError(
                        f"PSF aperture flux mismatch {step}_{nstar} {los} {band}: "
                        f"{aperture_from_particles} vs {aperture_from_image}")
                signed = signed_distance_arcsec(region)
                dist = np.full(len(flux), np.inf, dtype=np.float64)
                dist[on] = signed[iy[on], ix[on]]
                boundary_masks = {
                    width: on & (np.abs(dist) <= width + 1.0e-10)
                    for width in BOUNDARY_DISTANCE_THRESHOLDS_AS
                }
            else:
                q = np.zeros(len(flux), dtype=np.float64)
                dist = np.full(len(flux), np.inf, dtype=np.float64)
                boundary_masks = {width: np.zeros(len(flux), bool)
                                  for width in BOUNDARY_DISTANCE_THRESHOLDS_AS}

            w = flux[:, iband]
            primary = primary_b[(step, nstar, los, band)]
            pos_scores = {}
            psf_scores = {}
            boundary_scores = {width: {} for width in BOUNDARY_DISTANCE_THRESHOLDS_AS}
            for b in B_FINE:
                b = float(b)
                keep = keeps[b]
                pos = score_binary(w, keep, inside)
                psf = score_fractional(w, keep, q)
                pos_scores[b] = pos
                psf_scores[b] = psf
                for width in BOUNDARY_DISTANCE_THRESHOLDS_AS:
                    boundary_scores[width][b] = score_boundary(
                        w, keep, inside, boundary_masks[width])

                # Independently verify all nine published corrected-distance scores.
                if round(b, 2) in np.round(np.asarray(sb.B_GRID, float), 2):
                    ref = corrected_ref[(step, nstar, los, band, round(b, 2))]
                    for name, got in zip(("iou", "completeness", "precision"), pos):
                        if not (np.isnan(got) and np.isnan(ref[name])) and abs(got - ref[name]) > 1.2e-4:
                            raise RuntimeError(
                                f"Corrected positional {name} mismatch {step}_{nstar} "
                                f"{los} {band} b={b:.2f}: {got} != {ref[name]}")

            pos_best = float(B_FINE[int(np.nanargmax([pos_scores[float(b)][0] for b in B_FINE]))])
            psf_best = float(B_FINE[int(np.nanargmax([psf_scores[float(b)][0] for b in B_FINE]))])
            for b in B_FINE:
                b = float(b)
                pos = pos_scores[b]
                psf = psf_scores[b]
                row = {
                    "step": step, "nstar": nstar, "z": z, "los": los, "band": band,
                    "b": b, "photometry_factor": factor, "area_arcsec2": area,
                    "pos_iou": pos[0], "pos_completeness": pos[1], "pos_precision": pos[2],
                    "psf_iou": psf[0], "psf_completeness": psf[1], "psf_precision": psf[2],
                    "pos_fine_best_b": pos_best, "psf_fine_best_b": psf_best,
                    "primary_coarse_best_b": primary,
                }
                # Preserve the saved ``w`` column names for cache compatibility;
                # these values are pixel-centre distance thresholds.
                for width in BOUNDARY_DISTANCE_THRESHOLDS_AS:
                    score = boundary_scores[width][b]
                    token = str(width).replace(".", "p")
                    row[f"boundary_iou_w{token}_arcsec"] = score[0]
                    row[f"boundary_completeness_w{token}_arcsec"] = score[1]
                    row[f"boundary_precision_w{token}_arcsec"] = score[2]
                rows.append(row)
            print(f"    {los}/{band}: area={area:.3f} arcsec2; "
                  f"pos b*={pos_best:.2f}, PSF b*={psf_best:.2f}", flush=True)
            del image, labels, mask

    del raw, keeps, flux
    return rows


def quantiles(values: list[float]) -> tuple[float, float, float]:
    a = np.asarray(values, dtype=float)
    a = a[np.isfinite(a)]
    if a.size == 0:
        return np.nan, np.nan, np.nan
    return tuple(float(x) for x in np.quantile(a, (0.16, 0.50, 0.84)))


def make_summaries(rows: list[dict], n_halo: int) -> tuple[list[dict], list[dict], list[dict]]:
    by_halo_view_b = {(r["step"], int(r["nstar"]), r["los"], r["band"], round(float(r["b"]), 2)): r
                      for r in rows}
    bands = sorted({r["band"] for r in rows})
    parameter_rows = []
    curve_rows = []
    for band in bands:
        per_halo = defaultdict(list)
        per_b_metric = defaultdict(list)
        for r in rows:
            if r["band"] != band:
                continue
            key = (r["step"], int(r["nstar"]), round(float(r["b"]), 2))
            per_halo[key].append(r)
        halo_metrics = {"pos": {}, "psf": {}}
        for (step, nstar, b), views in per_halo.items():
            for metric, keyname in (("pos", "pos_iou"), ("psf", "psf_iou")):
                vals = [float(v[keyname]) for v in views]
                halo_metrics[metric].setdefault(b, []).append(float(np.nanmedian(vals)))
        for metric in ("pos", "psf"):
            selected = defaultdict(list)
            max_scores = defaultdict(list)
            for step, nstar in sorted({(r["step"], int(r["nstar"])) for r in rows
                                       if r["band"] == band}):
                view_best = []
                view_peak = []
                for los in rv.LOS:
                    view_rows = [r for r in rows if r["band"] == band and r["step"] == step
                                 and int(r["nstar"]) == nstar and r["los"] == los]
                    keyname = f"{metric}_iou"
                    best = max(view_rows, key=lambda r: (-np.inf if not np.isfinite(float(r[keyname]))
                                                         else float(r[keyname])))
                    view_best.append(float(best["b"]))
                    view_peak.append(float(best[keyname]))
                selected["b"].append(float(np.median(view_best)))
                max_scores["peak"].append(float(np.median(view_peak)))
            q16, med, q84 = quantiles(selected["b"])
            p16, pmed, p84 = quantiles(max_scores["peak"])
            score_key = f"{metric}_iou"
            at_010 = quantiles(halo_metrics[metric].get(0.10, []))[1]
            at_020 = quantiles(halo_metrics[metric].get(0.20, []))[1]
            parameter_rows.append({
                "band": band,
                "criterion": "positional_weighted_iou" if metric == "pos" else "psf_aperture_weighted_iou",
                "n_halo_outputs": n_halo,
                "median_host_b_across_3_view_medians": med,
                "p16_host_b": q16, "p84_host_b": q84,
                "median_host_view_optimised_iou": pmed,
                "p16_host_optimised_iou": p16, "p84_host_optimised_iou": p84,
                "median_iou_at_b0.10": at_010, "median_iou_at_b0.20": at_020,
            })
            for b, vals in sorted(halo_metrics[metric].items()):
                lo, m, hi = quantiles(vals)
                curve_rows.append({"band": band, "criterion": metric,
                                   "b": b, "n_halo_outputs": len(vals),
                                   "p16_host_median_iou": lo,
                                   "median_host_median_iou": m,
                                   "p84_host_median_iou": hi})

    boundary_rows = []
    for band in bands:
        for width in BOUNDARY_DISTANCE_THRESHOLDS_AS:
            token = str(width).replace(".", "p")
            for choice in ("primary_coarse_best", "b0.20"):
                host_vals = []
                for step, nstar in sorted({(r["step"], int(r["nstar"])) for r in rows
                                           if r["band"] == band}):
                    view_vals = []
                    for los in rv.LOS:
                        if choice == "b0.20":
                            b = 0.20
                        else:
                            b = float(by_halo_view_b[(step, nstar, los, band, 0.10)]
                                      ["primary_coarse_best_b"])
                        record = by_halo_view_b[(step, nstar, los, band, round(b, 2))]
                        view_vals.append(float(record[f"boundary_iou_w{token}_arcsec"]))
                    host_vals.append(float(np.nanmedian(view_vals)))
                lo, med, hi = quantiles(host_vals)
                boundary_rows.append({"band": band, "distance_threshold_arcsec": width,
                                      "criterion": choice, "n_halo_outputs": len(host_vals),
                                      "p16_host_median_iou": lo, "median_host_median_iou": med,
                                      "p84_host_median_iou": hi})
    # Add the calibrated-minus-0.20 difference at each boundary width.
    lookup = {(r["band"], r["distance_threshold_arcsec"], r["criterion"]): r
              for r in boundary_rows}
    for band in bands:
        for width in BOUNDARY_DISTANCE_THRESHOLDS_AS:
            a = lookup[(band, width, "primary_coarse_best")]
            b = lookup[(band, width, "b0.20")]
            a["median_difference_calibrated_minus_b0.20"] = (
                float(a["median_host_median_iou"]) - float(b["median_host_median_iou"]))
    return parameter_rows, curve_rows, boundary_rows


def make_plot(curve_rows: list[dict], outpath: Path) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.rcParams.update({"xtick.direction": "in", "ytick.direction": "in"})
    fig, axes = plt.subplots(1, 2, figsize=(8.3, 3.35), sharey=True)
    colors = {"pos": "#2468a2", "psf": "#c5653e"}
    labels = {"pos": "Particle-position weights", "psf": "PSF-aperture weights"}
    for ax, band in zip(axes, ("V", "g")):
        for criterion in ("pos", "psf"):
            subset = sorted((r for r in curve_rows if r["band"] == band
                             and r["criterion"] == criterion), key=lambda r: float(r["b"]))
            x = [float(r["b"]) for r in subset]
            y = [float(r["median_host_median_iou"]) for r in subset]
            low = [float(r["p16_host_median_iou"]) for r in subset]
            high = [float(r["p84_host_median_iou"]) for r in subset]
            ax.plot(x, y, marker="o", ms=3.4, lw=1.35, color=colors[criterion],
                    label=labels[criterion])
            ax.fill_between(x, low, high, color=colors[criterion], alpha=0.14, linewidth=0)
        ax.text(0.04, 0.96, band, transform=ax.transAxes, ha="left", va="top",
                fontsize=9, fontweight="bold")
        ax.set_xlabel(r"Linking parameter $b_{\rm link}$")
        ax.set_xlim(0.015, 0.205)
        ax.set_xticks(np.arange(0.025, 0.201, 0.025))
        ax.grid(alpha=0.2, linewidth=0.6)
        ax.tick_params(top=True, bottom=True, left=True, right=True, direction="in")
    axes[0].set_ylabel("Median light-weighted overlap")
    axes[1].legend(loc="lower right", frameon=False, fontsize=9)
    # Leave the seam's tick mark visible while suppressing only the first
    # numeric label on the right panel, which would collide with 0.200.
    axes[1].set_xticks(np.arange(0.025, 0.201, 0.025),
                       ["", *[f"{x:.3f}" for x in np.arange(0.050, 0.201, 0.025)]])
    fig.subplots_adjust(left=.12, right=.99, bottom=.18, top=.97, wspace=0.0)
    fig.savefig(outpath)
    plt.close(fig)


def make_boundary_plot(boundary_rows: list[dict], outpath: Path) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.rcParams.update({"xtick.direction": "in", "ytick.direction": "in"})
    colors = {"primary_coarse_best": "#2468a2", "b0.20": "#c5653e"}
    labels = {"primary_coarse_best": "Position-calibrated $b_{\\rm link}$",
              "b0.20": "$b_{\\rm link}=0.20$"}
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.15), sharey=True)
    for ax, band in zip(axes, ("V", "g")):
        for criterion in ("primary_coarse_best", "b0.20"):
            subset = sorted((r for r in boundary_rows
                             if r["band"] == band and r["criterion"] == criterion),
                            key=lambda r: float(r["distance_threshold_arcsec"]))
            x = np.asarray([float(r["distance_threshold_arcsec"]) for r in subset])
            med = np.asarray([float(r["median_host_median_iou"]) for r in subset])
            lo = np.asarray([float(r["p16_host_median_iou"]) for r in subset])
            hi = np.asarray([float(r["p84_host_median_iou"]) for r in subset])
            ax.plot(x, med, marker="o", ms=3.8, lw=1.4,
                    color=colors[criterion], label=labels[criterion])
            ax.fill_between(x, lo, hi, color=colors[criterion], alpha=.13, linewidth=0)
        ax.axvline(.5, color="0.3", lw=.8, ls="--")
        ax.set_xscale("log")
        ax.set_xticks([.2, .3, .5, 1., 2.], ["0.2", "0.3", "0.5", "1", "2"])
        ax.set_xlabel(r"Distance threshold $\Delta$ (arcsec)")
        ax.set_ylim(0., 1.)
        ax.grid(alpha=.2, linewidth=.6)
        ax.text(.04, .96, band, transform=ax.transAxes, ha="left", va="top",
                fontsize=9, fontweight="bold",
                bbox=dict(facecolor="white", edgecolor="none", alpha=.85, pad=1.))
        ax.tick_params(top=True, bottom=True, left=True, right=True, direction="in")
    axes[0].set_ylabel("Median overlap")
    axes[1].legend(loc="lower right", frameon=False, fontsize=9)
    # Keep the tick at the shared edge but avoid colliding with the 2-arcsec
    # label at the right edge of the V panel.
    axes[1].set_xticks([0.2, 0.3, 0.5, 1.0, 2.0], ["", "0.3", "0.5", "1", "2"])
    fig.subplots_adjust(left=.13, right=.99, bottom=.19, top=.97, wspace=0.0)
    fig.savefig(outpath)
    plt.close(fig)


def process_sample(sample: dict[str, str], corrected_ref: dict, primary_b: dict,
                   factors: dict, outdir: Path, index: int, total: int) -> tuple[Path, list[dict]]:
    step, nstar = sample["step"], int(sample["nstar"])
    print(f"[{index}/{total}] corrected PSF/boundary analysis {step}_{nstar}", flush=True)
    rows = measure_halo(step, nstar, corrected_ref, primary_b, factors)
    path = outdir / "per_halo" / f"{step}_{nstar}.tsv"
    write_tsv(path, rows)
    return path, rows


def write_summary_outputs(outdir: Path, rows: list[dict], n_halo: int) -> None:
    parameter_rows, curve_rows, boundary_rows = make_summaries(rows, n_halo)
    write_tsv(outdir / "fine_grid_parameter_summary.tsv", parameter_rows)
    write_tsv(outdir / "fine_grid_curve_summary.tsv", curve_rows)
    write_tsv(outdir / "boundary_width_summary.tsv", boundary_rows)
    make_plot(curve_rows, outdir / "corrected_positional_and_psf_overlap.pdf")
    make_boundary_plot(boundary_rows, outdir / "boundary_width_sensitivity.pdf")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=RESULTS)
    parser.add_argument("--limit", type=int, default=None,
                        help="Process only the first N sample entries (for a local preflight).")
    parser.add_argument("--only", help="Single sample entry STEP:NSTAR.")
    parser.add_argument("--workers", type=int, default=8,
                        help="Independent halo calculations to run concurrently (default: 8).")
    parser.add_argument("--resume", action="store_true",
                        help="skip completed per-halo files in an existing output directory")
    parser.add_argument("--summaries-only", action="store_true",
                        help="Rebuild summary tables and plot from an existing per-view results file.")
    args = parser.parse_args()
    outdir = args.output.resolve()
    if args.summaries_only:
        if not (outdir / "corrected_scores_per_view.tsv").is_file():
            raise FileNotFoundError(outdir / "corrected_scores_per_view.tsv")
        all_rows = read_tsv(outdir / "corrected_scores_per_view.tsv")
        n_halo = sum(1 for r in read_tsv(MANIFEST)
                     if r["calibration_retained"].lower() == "true")
        write_summary_outputs(outdir, all_rows, n_halo)
        print(f"Rebuilt summaries from {len(all_rows)} saved rows in {outdir}", flush=True)
        return
    if outdir.exists() and any(outdir.iterdir()) and not args.resume:
        raise FileExistsError(f"Refusing to overwrite non-empty analysis directory: {outdir}")
    outdir.mkdir(parents=True, exist_ok=True)

    corrected_ref, primary_b, factors = load_references()
    samples = [r for r in read_tsv(MANIFEST) if r["calibration_retained"].lower() == "true"]
    if args.only:
        only_step, only_nstar = args.only.split(":", 1)
        samples = [r for r in samples if r["step"] == only_step.zfill(5)
                   and int(r["nstar"]) == int(only_nstar)]
        if len(samples) != 1:
            raise RuntimeError(f"Expected one sample for {args.only}; found {len(samples)}")
    if args.limit is not None:
        samples = samples[:args.limit]
    if not samples:
        raise RuntimeError("No calibration sample entries selected")

    all_rows: list[dict] = []
    per_halo_dir = outdir / "per_halo"
    per_halo_dir.mkdir(exist_ok=True)
    pending = [sample for sample in samples
               if not (per_halo_dir / f"{sample['step']}_{int(sample['nstar'])}.tsv").is_file()]
    if args.resume:
        print(f"Resuming: {len(samples) - len(pending)} of {len(samples)} halo entries already saved; "
              f"{len(pending)} remain", flush=True)
    if args.workers <= 1 or len(pending) == 1:
        for i, sample in enumerate(pending, start=1):
            _, rows = process_sample(sample, corrected_ref, primary_b, factors,
                                     outdir, i, len(pending))
            all_rows.extend(rows)
    elif pending:
        context = mp.get_context("fork")
        with ProcessPoolExecutor(max_workers=args.workers, mp_context=context) as pool:
            futures = {
                pool.submit(process_sample, sample, corrected_ref, primary_b, factors,
                            outdir, i, len(pending)): i
                for i, sample in enumerate(pending, start=1)
            }
            finished = {}
            for future in as_completed(futures):
                index = futures[future]
                path, rows = future.result()
                finished[index] = rows
                print(f"[{index}/{len(pending)}] saved {path.name}", flush=True)
            for i in sorted(finished):
                all_rows.extend(finished[i])

    # Re-read every per-halo file so a resumed run produces the same complete
    # summary as a single uninterrupted run.
    all_rows = []
    for sample in samples:
        path = per_halo_dir / f"{sample['step']}_{int(sample['nstar'])}.tsv"
        if not path.is_file():
            raise FileNotFoundError(f"Missing completed halo result {path}")
        all_rows.extend(read_tsv(path))

    write_tsv(outdir / "corrected_scores_per_view.tsv", all_rows)
    write_summary_outputs(outdir, all_rows, len(samples))
    manifest = {
        "sample_entries": len(samples),
        "photometry": "corrected luminosity-distance factor applied before thresholding",
        "grid": B_FINE.tolist(),
        "boundary_distance_thresholds_arcsec": list(BOUNDARY_DISTANCE_THRESHOLDS_AS),
        "boundary_definition": "absolute signed Euclidean distance from each pixel centre to the nearest pixel centre of the opposite isophotal class; pixel scale 0.2 arcsec",
        "legacy_one_pixel_ring_equivalent_distance_threshold_arcsec": 0.3,
        "psf_metric": "each star receives the fraction of its normalized PSF kernel inside the corrected isophotal region",
        "coarse_membership_source": str(MEMBERSHIP_ROOT),
        "fine_nodes_computed": [0.06, 0.07, 0.09, 0.11, 0.13, 0.14],
        "corrected_positional_validation_rows": sum(
            1 for r in all_rows if round(float(r["b"]), 2) in np.round(np.asarray(sb.B_GRID, float), 2)),
    }
    (outdir / "analysis_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(f"WROTE {len(all_rows)} per-view/group rows to {outdir}", flush=True)


if __name__ == "__main__":
    main()
