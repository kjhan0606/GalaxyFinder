#!/usr/bin/env python3
"""Render and score the TNG300 SubhaloPos-centred calibration follow-up.

Each task handles one preselected FoF host and uses the catalogue position of
the maximally bound particle for the image and aperture centre. The
stellar-density peak remains the FoF root. No linking parameter is tuned on
TNG.
"""
from __future__ import annotations

import argparse
import csv
import fcntl
import os
import sys
from pathlib import Path

os.environ.setdefault("SPS_HOME", "/home/kjhan/fsps_data")
for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[name] = "1"

import h5py
import numpy as np
from astropy.cosmology import FlatLambdaCDM
from scipy.signal import fftconvolve

HERE = Path(__file__).resolve().parent
SOURCE = HERE.parent.parent / "paper_prepublication_candidate_20261028/supplement/tng_transfer_protocol/source_dependencies"
sys.path.insert(0, str(SOURCE))
import sb_vs_b as sb  # noqa: E402
import sb_revise as rv  # noqa: E402

OUT = HERE / "rendered_scores"
DATA = HERE.parent / "tng300_blind_transfer_20261029/host_star_membership"
SAMPLE = HERE.parent / "tng300_blind_transfer_20261029/tng300_massmatched_40_hosts.csv"
CENTERS = HERE / "tng300_subhalo_pos_40.csv"
Z_SNAP = 0.19728418237600986
A_SNAP = 1.0 / (1.0 + Z_SNAP)
H_TNG = 0.6774
OM_TNG = 0.3089
BOX_CKPC_H = 205000.0
COSMO_TNG = FlatLambdaCDM(H0=67.74, Om0=OM_TNG)
COSMO_FSPS = FlatLambdaCDM(H0=72.0, Om0=0.27)
METHODS = ("b020", "b010", "cal_primary", "cal_lowz_loo", "aperture30")
BAND_INDEX = {"V": 0, "g": 1}


def read_sample(rank: int) -> dict[str, str]:
    with SAMPLE.open(newline="") as f:
        rows = [r for r in csv.DictReader(f) if int(r["sample_rank"]) == rank]
    if len(rows) != 1:
        raise RuntimeError(f"Expected one locked sample row for rank {rank}; got {len(rows)}")
    return rows[0]


def time_grids(a_snap: float):
    """RAMSES-style birth-time mapping evaluated with the TNG cosmology."""
    a = np.linspace(0.02, a_snap, 4000)
    grid = np.linspace(a_snap, 1.0, 4000)
    e = lambda q: np.sqrt(OM_TNG / q**3 + (1.0 - OM_TNG))
    tau_max = -np.trapezoid(1.0 / (grid**3 * e(grid)), grid)
    inv_tau = 1.0 / (a**3 * e(a))
    integ = np.zeros_like(a)
    integ[1:] = np.cumsum(0.5 * (inv_tau[:-1] + inv_tau[1:]) * np.diff(a))
    tau = (tau_max - integ[-1]) + integ
    inv_t = 1.0 / (a * e(a))
    tcum = np.zeros_like(a)
    tcum[1:] = np.cumsum(0.5 * (inv_t[:-1] + inv_t[1:]) * np.diff(a)) * (977.792221 / 67.74)
    return a, tau, tcum


def distance_factor(z: float) -> float:
    """Convert FSPS's reference distance to the TNG cosmological distance."""
    dl_fsps = float(COSMO_FSPS.luminosity_distance(z).value)
    dl_tng = float(COSMO_TNG.luminosity_distance(z).value)
    return (dl_fsps / dl_tng) ** 2


def get_fsps_grid():
    cache = OUT / "fsps_cache" / f"fsps_tng_z{Z_SNAP:.6f}.npz"
    cache.parent.mkdir(parents=True, exist_ok=True)
    with (OUT / "fsps_cache" / "grid.lock").open("w") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        if cache.exists():
            with np.load(cache) as d:
                if abs(float(d["z"]) - Z_SNAP) > 1e-12:
                    raise RuntimeError("Cached FSPS grid has a different redshift")
                grid = (d["age"].copy(), d["logz"].copy(), d["flux"].copy())
        else:
            grid = sb.ssp_flux(Z_SNAP, cache)
        fcntl.flock(lock, fcntl.LOCK_UN)
    return grid


def wrapped_center(coords: np.ndarray, box: float) -> np.ndarray:
    """Periodic circular mean, used only to place the host in the image stamp."""
    angle = coords * (2.0 * np.pi / box)
    return (np.arctan2(np.sin(angle).mean(axis=0), np.cos(angle).mean(axis=0)) % (2.0 * np.pi)) * box / (2.0 * np.pi)


def read_pos(group_id: int, central_id: int) -> np.ndarray:
    if not CENTERS.is_file():
        raise FileNotFoundError(f"Missing verified SubhaloPos manifest: {CENTERS}")
    with CENTERS.open(newline="") as f:
        matches = [r for r in csv.DictReader(f) if int(r["group_id"]) == group_id]
    if len(matches) != 1 or int(matches[0]["central_subhalo_id"]) != central_id:
        raise RuntimeError(f"SubhaloPos manifest does not match group {group_id} central {central_id}")
    row = matches[0]
    return np.array([float(row[f"pos_{axis}_ckpc_h"]) for axis in "xyz"], dtype=np.float64)


def load_host(rank: int):
    sample = read_sample(rank)
    gid = int(sample["group_id"])
    path = DATA / f"halo_{gid:08d}_stars_membership.h5"
    if not path.is_file():
        raise FileNotFoundError(path)
    with h5py.File(path, "r") as f:
        g = f["PartType4"]
        coord_ckpc_h = g["Coordinates"][:].astype(np.float64)
        mass_tng = g["Masses"][:].astype(np.float64)
        initial_tng = g["GFM_InitialMass"][:].astype(np.float64)
        metallicity = g["GFM_Metallicity"][:].astype(np.float64)
        birth_a = g["GFM_StellarFormationTime"][:].astype(np.float64)
        removed_satellite = g["RemoveSatelliteGalaxy"][:].astype(bool)
        z = float(f.attrs["Redshift"])
        box_ckpc_h = float(f.attrs.get("BoxSize", BOX_CKPC_H))
        central_id = int(f.attrs["CentralSubhaloID"])
    if abs(z - Z_SNAP) > 1.0e-6:
        raise RuntimeError(f"FoF group {gid}: snapshot redshift {z} differs from locked {Z_SNAP}")

    # TNG PartType4 also contains wind-phase gas cells (formation scale factor
    # <= 0); those are not stars and must not enter images or stellar FoF.
    valid_star = np.isfinite(birth_a) & (birth_a > 0.0) & (birth_a <= 1.0 / (1.0 + z))
    keep_parent = (~removed_satellite) & valid_star
    if keep_parent.sum() < 2:
        raise RuntimeError(f"FoF group {gid}: fewer than two retained real stars")
    xyz = coord_ckpc_h[keep_parent] / 1000.0  # comoving Mpc/h
    center = read_pos(gid, central_id) / 1000.0  # SubhaloPos, comoving Mpc/h
    delta = xyz - center
    delta -= (box_ckpc_h / 1000.0) * np.round(delta / (box_ckpc_h / 1000.0))

    # Keep the density-peak star as the FoF root, with coordinates relative
    # to the mass-weighted catalogue centre for images and apertures.
    seed = sb.nucleus_index(delta[:, 0], delta[:, 1], delta[:, 2])
    seed_offset_pkpc = (
        np.linalg.norm(delta[seed]) * (1.0 / (1.0 + z)) / H_TNG * 1000.0
    )
    current_mass_msun_h = mass_tng[keep_parent] * 1.0e10
    initial_mass_msun = initial_tng[keep_parent] * 1.0e10 / H_TNG
    a_grid, tau_grid, t_cum = time_grids(1.0 / (1.0 + z))
    birth_a = np.clip(birth_a[keep_parent], a_grid[0], 1.0 / (1.0 + z))
    birth_tau = np.interp(birth_a, a_grid, tau_grid)
    return (sample, gid, z, delta, current_mass_msun_h, initial_mass_msun,
            birth_tau, metallicity[keep_parent], valid_star, removed_satellite,
            seed, seed_offset_pkpc, center * 1000.0, (a_grid, tau_grid, t_cum))


def method_masks(dx, dy, dz, mass, center_idx):
    # The density-peak seed defines the linked component; positions remain
    # centred on SubhaloPos for the image and aperture measurements.
    masks = {}
    for b in (0.09, 0.10, 0.12, 0.20):
        link = b * np.cbrt(mass / (2.7755e11 * OM_TNG))
        root = sb.fof_largest(dx, dy, dz, link)
        masks[b] = root == root[center_idx]
    return masks


def weighted_r50(radius_kpc, weights, selection):
    """Return the projected radius enclosing half of the selected full flux."""
    use = selection & np.isfinite(weights) & (weights > 0.0)
    if not np.any(use):
        return float("nan")
    radius = radius_kpc[use]
    weight = weights[use]
    order = np.argsort(radius)
    cumulative = np.cumsum(weight[order])
    return float(radius[order[np.searchsorted(cumulative, 0.5 * cumulative[-1], side="left")]])


def render_one(rank: int):
    OUT.mkdir(parents=True, exist_ok=True)
    outpath = OUT / f"host_rank_{rank:02d}.csv"
    if outpath.exists():
        raise FileExistsError(f"Refusing to overwrite existing score file: {outpath}")
    (sample, gid, z, xyz, mass, initial_mass, birth_tau, metallicity,
     valid_star, removed_satellite, seed, seed_offset_pkpc, center_ckpc_h,
     time_grid) = load_host(rank)
    print(f"START rank={rank} FoF={gid} retained={xyz.shape[0]} "
          f"z={z:.6f} seed_offset_pkpc={seed_offset_pkpc:.3f}", flush=True)
    age_grid, logz_grid, ssp = get_fsps_grid()
    a_grid, tau_grid, t_cum = time_grid
    flux = sb.star_flux(birth_tau, metallicity, initial_mass, 1.0 / (1.0 + z),
                        a_grid, tau_grid, t_cum, age_grid, logz_grid, ssp)
    flux *= distance_factor(z)
    masks = method_masks(xyz[:, 0], xyz[:, 1], xyz[:, 2], mass, seed)
    dm = float(COSMO_TNG.comoving_transverse_distance(z).value)
    scale = 206265.0 / (H_TNG * dm)  # arcsec per comoving Mpc/h
    n = int(round(2.0 * sb.STAMP_AS / sb.PIX_AS))
    kernels = [sb.psf_kernel(sb.PIX_AS, sb.FWHM_AS, sb.MOFFAT_BETA, wing, sb.KERNEL_AS)[0]
               for wing in sb.WING_A]
    rows = []
    axes = ((1, 2, "x"), (2, 0, "y"), (0, 1, "z"))
    calibration_b = {
        "cal_primary": {"V": 0.12, "g": 0.10},
        "cal_lowz_loo": {"V": 0.10, "g": 0.09},
    }
    for i_view, (ia, ib, los) in enumerate(axes):
        radius_kpc = np.hypot(xyz[:, ia], xyz[:, ib]) * (1.0 / (1.0 + z)) / H_TNG * 1000.0
        ix = np.rint(xyz[:, ia] * scale / sb.PIX_AS).astype(np.int32) + n // 2
        iy = np.rint(xyz[:, ib] * scale / sb.PIX_AS).astype(np.int32) + n // 2
        on = (ix >= 0) & (ix < n) & (iy >= 0) & (iy < n)
        masks_view = dict(masks)
        masks_view["aperture30"] = (
            np.hypot(xyz[:, ia], xyz[:, ib]) / H_TNG
            * (1.0 / (1.0 + z)) * 1000.0 <= 30.0)
        for iband, band in enumerate(("V", "g")):
            observed_raw = sb.bin_flux(ix, iy, flux[:, iband], n)
            observed_img = fftconvolve(observed_raw, kernels[iband], mode="same")
            observed_img[observed_img < 0.0] = 0.0
            observed_region, observed_area = sb.central_isophote(
                observed_img, sb.THRESHOLDS[iband], sb.PIX_AS)
            inside = np.zeros(ix.size, bool)
            inside[on] = observed_region[iy[on], ix[on]]
            observed_edge = bool(observed_region[0, :].any() or observed_region[-1, :].any()
                                 or observed_region[:, 0].any() or observed_region[:, -1].any())
            l_observed = float(flux[inside, iband].sum())
            r_observed = float(np.sqrt(observed_area / np.pi)) if observed_area > 0 else np.nan
            r50_observed = weighted_r50(radius_kpc, flux[:, iband], inside)
            method_bs = {
                "b020": 0.20,
                "b010": 0.10,
                "cal_primary": calibration_b["cal_primary"][band],
                "cal_lowz_loo": calibration_b["cal_lowz_loo"][band],
            }
            for method in METHODS:
                if method == "aperture30":
                    keep = masks_view[method]
                else:
                    keep = masks_view[method_bs[method]]
                selected_raw = sb.bin_flux(ix, iy, flux[:, iband] * keep, n)
                selected_img = fftconvolve(selected_raw, kernels[iband], mode="same")
                selected_img[selected_img < 0.0] = 0.0
                selected_region, selected_area = sb.central_isophote(
                    selected_img, sb.THRESHOLDS[iband], sb.PIX_AS)
                selected_edge = bool(selected_region[0, :].any() or selected_region[-1, :].any()
                                     or selected_region[:, 0].any() or selected_region[:, -1].any())
                iou, completeness, precision = rv.light_scores(flux[:, iband], keep, inside)
                l_selected = float(flux[keep, iband].sum())
                r_selected = float(np.sqrt(selected_area / np.pi)) if selected_area > 0 else np.nan
                r50_selected = weighted_r50(radius_kpc, flux[:, iband], keep)
                log_l = (float(np.log10(l_selected / l_observed))
                         if l_selected > 0 and l_observed > 0 else np.nan)
                log_r = (float(np.log10(r_selected / r_observed))
                         if r_selected > 0 and r_observed > 0 else np.nan)
                log_r50 = (float(np.log10(r50_selected / r50_observed))
                           if r50_selected > 0 and r50_observed > 0 else np.nan)
                contamination = 1.0 - float(precision) if np.isfinite(precision) else np.nan
                rows.append({
                    "sample_rank": rank, "group_id": gid, "redshift": z,
                    "los": los, "band": band, "method": method,
                    "b": method_bs.get(method, np.nan),
                    "image_center_rule": "TNG SubhaloPos of selected central subhalo",
                    "image_center_x_ckpc_h": float(center_ckpc_h[0]),
                    "image_center_y_ckpc_h": float(center_ckpc_h[1]),
                    "image_center_z_ckpc_h": float(center_ckpc_h[2]),
                    "density_seed_offset_pkpc": seed_offset_pkpc,
                    "psf_fwhm_arcsec": float(sb.FWHM_AS),
                    "n_parent_type4": int(valid_star.size),
                    "n_retained_real_stars": int(xyz.shape[0]),
                    "n_nonpositive_birth_time": int((~valid_star).sum()),
                    "n_removed_satellite_type4": int(removed_satellite.sum()),
                    "observed_area_arcsec2": observed_area,
                    "selected_area_arcsec2": selected_area,
                    "observed_edge": int(observed_edge),
                    "observed_missing": int(observed_area <= 0.0),
                    "selected_edge": int(selected_edge),
                    "isophote_missing": int(selected_area <= 0.0),
                    "iou": float(iou) if np.isfinite(iou) else np.nan,
                    "completeness": float(completeness) if np.isfinite(completeness) else np.nan,
                    "precision": float(precision) if np.isfinite(precision) else np.nan,
                    "contamination": contamination,
                    "log10_luminosity_ratio": log_l,
                    "log10_radius_ratio": log_r,
                    "log10_R50_ratio": log_r50,
                    "L_selected": l_selected, "L_observed": l_observed,
                    "R_selected_arcsec": r_selected, "R_observed_arcsec": r_observed,
                    "R50_selected_pkpc": r50_selected,
                    "R50_observed_pkpc": r50_observed,
                })
        print(f"  projected view {los} done", flush=True)

    fields = list(rows[0])
    with outpath.open("x", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    print(f"DONE rank={rank} FoF={gid} rows={len(rows)} output={outpath}", flush=True)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--sample-rank", type=int, required=True)
    args = p.parse_args()
    if not 1 <= args.sample_rank <= 40:
        raise SystemExit("sample rank must be in [1, 40]")
    render_one(args.sample_rank)


if __name__ == "__main__":
    main()
