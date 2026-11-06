#!/usr/bin/env python3
"""Measure stellar masses and projected half-light radii across TNG trees.

For each fixed SubLink progenitor state, the script renders mock images at
three observer redshifts. This cross-rendering holds the particle structure
fixed while changing distance, redshifted spectra, and cosmological dimming.
It compares an isophotal particle selection with fixed 3D stellar groups,
projected apertures, and Darwin-calibrated spherical apertures.
"""
from __future__ import annotations

import argparse
import csv
import fcntl
import os
import sys
from pathlib import Path

for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[name] = "1"

import h5py
import numpy as np
from astropy.cosmology import FlatLambdaCDM
from scipy.signal import fftconvolve

HERE = Path(__file__).resolve().parent
SOURCE = HERE / "source_dependencies"
sys.path.insert(0, str(SOURCE))
import sb_vs_b as sb  # noqa: E402
import sb_revise as rv  # noqa: E402

OUT = HERE / "growth_scores"
H_TNG = 0.6774
OM_TNG = 0.3089
BOX_CKPC_H = 205000.0
COSMO_TNG = FlatLambdaCDM(H0=67.74, Om0=OM_TNG)
COSMO_FSPS = FlatLambdaCDM(H0=72.0, Om0=0.27)
TARGET_SNAPS = (84, 67, 50)
OBSERVER_REDSHIFTS = (0.19728418237601, 0.503047523244883, 0.99729422578194)
SPHERICAL_RADIUS_KPC = {"V": 80.0, "g": 60.0}
AXES = ((1, 2, "x"), (2, 0, "y"), (0, 1, "z"))


def read_rows(path: Path):
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream))


def locked_cohort():
    return {int(row["sample_rank"]) for row in read_rows(HERE / "lineage_cohort.csv")
            if int(row["central_at_all_targets"]) == 1}


def epoch_rows(rank: int):
    rows = [row for row in read_rows(HERE / "lineage_epochs.csv")
            if int(row["sample_rank"]) == rank and int(row["snapshot"]) in TARGET_SNAPS]
    result = {int(row["snapshot"]): row for row in rows}
    if set(result) != set(TARGET_SNAPS):
        raise RuntimeError(f"Rank={rank} is missing a target SubLink epoch")
    if any(int(row["is_fof_central"]) != 1 for row in result.values()):
        raise RuntimeError(f"Rank={rank} is not central at every target epoch")
    return result


def membership_path(snapshot: int, rank: int, row: dict[str, str]) -> Path:
    if snapshot == 84:
        group_id = int(row["descendant_group_id"])
        path = HERE / "tree_epoch_membership" / (
            f"snap084_rank{rank:02d}_group{group_id:06d}.h5")
    else:
        group_id = int(row["progenitor_fof_group_id"])
        path = HERE / "tree_epoch_membership" / (
            f"snap{snapshot:03d}_rank{rank:02d}_group{group_id:06d}.h5")
    return path


def time_grids(a_snap: float):
    a = np.linspace(0.02, a_snap, 4000)
    grid = np.linspace(a_snap, 1.0, 4000)
    expansion = lambda q: np.sqrt(OM_TNG / q**3 + (1.0 - OM_TNG))
    tau_max = -np.trapezoid(1.0 / (grid**3 * expansion(grid)), grid)
    inverse_tau = 1.0 / (a**3 * expansion(a))
    integ = np.zeros_like(a)
    integ[1:] = np.cumsum(0.5 * (inverse_tau[:-1] + inverse_tau[1:]) * np.diff(a))
    tau = (tau_max - integ[-1]) + integ
    inverse_t = 1.0 / (a * expansion(a))
    tcum = np.zeros_like(a)
    tcum[1:] = np.cumsum(0.5 * (inverse_t[:-1] + inverse_t[1:]) * np.diff(a)) * (977.792221 / 67.74)
    return a, tau, tcum


def distance_factor(z: float) -> float:
    return (float(COSMO_FSPS.luminosity_distance(z).value)
            / float(COSMO_TNG.luminosity_distance(z).value)) ** 2


def get_fsps_grid(z: float):
    cache = OUT / "fsps_cache" / f"tng_ssp_observer_z{z:.6f}.npz"
    cache.parent.mkdir(parents=True, exist_ok=True)
    with (OUT / "fsps_cache" / "grid.lock").open("w") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        if cache.exists():
            with np.load(cache) as data:
                if abs(float(data["z"]) - z) > 1e-10:
                    raise RuntimeError(f"FSPS cache redshift mismatch at {z}")
                grid = data["age"].copy(), data["logz"].copy(), data["flux"].copy()
        else:
            grid = sb.ssp_flux(z, cache)
        fcntl.flock(lock, fcntl.LOCK_UN)
    return grid


def weighted_r50(radius_pkpc, weights, selection):
    use = selection & np.isfinite(weights) & (weights > 0.0)
    if not np.any(use):
        return float("nan")
    radius = radius_pkpc[use]
    weight = weights[use]
    order = np.argsort(radius)
    cumulative = np.cumsum(weight[order])
    return float(radius[order[np.searchsorted(cumulative, 0.5 * cumulative[-1], side="left")]])


def load_epoch(rank: int, snapshot: int, row: dict[str, str]):
    path = membership_path(snapshot, rank, row)
    if not path.is_file():
        raise FileNotFoundError(path)
    with h5py.File(path, "r") as handle:
        data = handle["PartType4"]
        coord = data["Coordinates"][:].astype(np.float64)
        current_mass = data["Masses"][:].astype(np.float64)
        initial_mass = data["GFM_InitialMass"][:].astype(np.float64)
        metallicity = data["GFM_Metallicity"][:].astype(np.float64)
        birth_a = data["GFM_StellarFormationTime"][:].astype(np.float64)
        remove_satellite = data["RemoveSatelliteGalaxy"][:].astype(bool)
        central_subhalo_particle = data["CentralSubhaloParticle"][:].astype(bool)
        z_phys = float(handle.attrs["Redshift"])
        box_ckpc_h = float(handle.attrs.get("BoxSize", BOX_CKPC_H))
        group_id = int(handle.attrs["GroupID"])
        central_id = int(handle.attrs["CentralSubhaloID"])
    if group_id != int(row["progenitor_fof_group_id"]):
        raise RuntimeError(f"FoF group mismatch for rank={rank}, snapshot={snapshot}")
    if central_id != int(row["progenitor_subhalo_id"]):
        raise RuntimeError(f"SubLink central mismatch for rank={rank}, snapshot={snapshot}")
    if central_subhalo_particle.size != coord.shape[0] or remove_satellite.size != coord.shape[0]:
        raise RuntimeError(f"Membership length mismatch for rank={rank}, snapshot={snapshot}")
    if np.any(central_subhalo_particle & remove_satellite):
        raise RuntimeError(f"Overlapping central and satellite masks for rank={rank}, snapshot={snapshot}")

    valid_star = (np.isfinite(birth_a) & (birth_a > 0.0)
                  & (birth_a <= 1.0 / (1.0 + z_phys)))
    use = (~remove_satellite) & valid_star
    if not np.any(use & central_subhalo_particle):
        raise RuntimeError(f"No valid bound central stars at rank={rank}, snapshot={snapshot}")
    center = np.array([float(row[f"subhalo_{axis}_ckpc_h"]) for axis in "xyz"], dtype=float)
    delta = coord[use] / 1000.0 - center / 1000.0
    box = box_ckpc_h / 1000.0
    delta -= box * np.round(delta / box)
    mass_h = current_mass[use] * 1.0e10
    mass_msun = mass_h / H_TNG
    formed_mass_msun = initial_mass[use] * 1.0e10 / H_TNG
    metallicity = metallicity[use]
    birth_a = birth_a[use]
    central = central_subhalo_particle[use]
    seed = sb.nucleus_index(delta[:, 0], delta[:, 1], delta[:, 2])
    a_snap = 1.0 / (1.0 + z_phys)
    a_grid, tau_grid, t_cum = time_grids(a_snap)
    birth_a = np.clip(birth_a, a_grid[0], a_snap)
    birth_tau = np.interp(birth_a, a_grid, tau_grid)
    r3d_pkpc = np.linalg.norm(delta, axis=1) * a_snap / H_TNG * 1000.0
    link_masks = {}
    for b in (0.10, 0.12, 0.20):
        link = b * np.cbrt(mass_h / (2.7755e11 * OM_TNG))
        labels = sb.fof_largest(delta[:, 0], delta[:, 1], delta[:, 2], link)
        link_masks[b] = labels == labels[seed]
    ref_mass_cat = float(row["stellar_mass_1e10_msun_h"]) * 1.0e10 / H_TNG
    ref_mass_particles = float(mass_msun[central].sum())
    catalog_mass_ratio = ref_mass_particles / ref_mass_cat if ref_mass_cat > 0 else float("nan")
    expected_n = int(row["stellar_particle_count"])
    central_n = int(central_subhalo_particle.sum())
    if central_n != expected_n:
        raise RuntimeError(
            f"Central particle count {central_n} != SubLink catalog count {expected_n} "
            f"for rank={rank}, snapshot={snapshot}"
        )
    return {
        "rank": rank, "snapshot": snapshot, "row": row, "z_phys": z_phys,
        "delta_comoving_mpc_h": delta, "radius3d_pkpc": r3d_pkpc,
        "mass_msun": mass_msun, "current_mass_1e10_msun_h": mass_h,
        "formed_mass_msun": formed_mass_msun, "metallicity": metallicity,
        "birth_tau": birth_tau, "central_mask": central, "link_masks": link_masks,
        "seed": seed, "a_snap": a_snap, "a_grid": a_grid, "tau_grid": tau_grid,
        "t_cum": t_cum, "catalog_mass_ratio": catalog_mass_ratio,
        "catalog_mass_msun": ref_mass_cat, "particle_central_mass_msun": ref_mass_particles,
        "n_particles": int(delta.shape[0]), "central_particle_count": central_n,
    }


def render_epoch(data, observer_redshift: float, grid_cache: dict):
    rank = data["rank"]
    snapshot = data["snapshot"]
    row = data["row"]
    z_phys = data["z_phys"]
    if observer_redshift not in grid_cache:
        grid_cache[observer_redshift] = get_fsps_grid(observer_redshift)
    age_g, logz_g, ssp = grid_cache[observer_redshift]
    flux = sb.star_flux(data["birth_tau"], data["metallicity"], data["formed_mass_msun"],
                        data["a_snap"], data["a_grid"], data["tau_grid"], data["t_cum"],
                        age_g, logz_g, ssp)
    flux *= distance_factor(observer_redshift)

    z = observer_redshift
    da_mpc = float(COSMO_TNG.angular_diameter_distance(z).value)
    if abs(observer_redshift - z_phys) <= 1.0e-10:
        # Match the established TNG renderer's comoving-distance expression
        # at the physical redshift, including identical pixel rounding.
        dm_mpc = float(COSMO_TNG.comoving_transverse_distance(z).value)
        scale_arcsec_per_comoving_mpc_h = 206265.0 / (H_TNG * dm_mpc)
    else:
        # Keep the physical particle coordinates fixed while changing only
        # the observer distance in the counterfactual renders.
        scale_arcsec_per_comoving_mpc_h = (
            206265.0 * data["a_snap"] / (H_TNG * da_mpc)
        )
    n = int(round(2.0 * sb.STAMP_AS / sb.PIX_AS))
    axes = AXES
    kernels = [sb.psf_kernel(sb.PIX_AS, sb.FWHM_AS, sb.MOFFAT_BETA, wing, sb.KERNEL_AS)[0]
               for wing in sb.WING_A]
    xyz_pkpc = data["delta_comoving_mpc_h"] * data["a_snap"] / H_TNG * 1000.0
    rows = []
    for ia, ib, los in axes:
        radius2_pkpc = np.hypot(xyz_pkpc[:, ia], xyz_pkpc[:, ib])
        ix = np.rint(data["delta_comoving_mpc_h"][:, ia]
                     * scale_arcsec_per_comoving_mpc_h / sb.PIX_AS).astype(np.int32) + n // 2
        iy = np.rint(data["delta_comoving_mpc_h"][:, ib]
                     * scale_arcsec_per_comoving_mpc_h / sb.PIX_AS).astype(np.int32) + n // 2
        on_image = (ix >= 0) & (ix < n) & (iy >= 0) & (iy < n)
        common = {
            "sample_rank": rank,
            "descendant_group_id": int(row["descendant_group_id"]),
            "descendant_subhalo_id": int(row["descendant_subhalo_id"]),
            "snapshot": snapshot,
            "physical_redshift": z_phys,
            "observer_redshift": observer_redshift,
            "line_of_sight": los,
            "n_input_stars": data["n_particles"],
            "n_central_bound_stars": data["central_particle_count"],
            "central_catalog_particle_mass_ratio": data["catalog_mass_ratio"],
        }
        for iband, band in enumerate(("V", "g")):
            observed_raw = sb.bin_flux(ix, iy, flux[:, iband], n)
            observed_img = fftconvolve(observed_raw, kernels[iband], mode="same")
            observed_img[observed_img < 0.0] = 0.0
            observed_region, observed_area = sb.central_isophote(
                observed_img, sb.THRESHOLDS[iband], sb.PIX_AS)
            cy = cx = observed_region.shape[0] // 2
            central_pixel_detected = bool(observed_region[cy, cx])
            region_detected = bool(observed_area > 0.0)
            if central_pixel_detected:
                isophote_status = "central_pixel_detected"
            elif region_detected:
                isophote_status = "nearest_bright_pixel_fallback"
            else:
                isophote_status = "no_bright_region"
            edge_contact = bool(
                region_detected and (
                    observed_region[0, :].any() or observed_region[-1, :].any()
                    or observed_region[:, 0].any() or observed_region[:, -1].any()))
            if region_detected:
                ry, rx = np.nonzero(observed_region)
                region_centroid_offset_pix = float(np.hypot(rx.mean() - cx, ry.mean() - cy))
                nearest_region_pixel_offset_pix = float(np.sqrt(
                    np.min((rx - cx) ** 2 + (ry - cy) ** 2)))
            else:
                region_centroid_offset_pix = float("nan")
                nearest_region_pixel_offset_pix = float("nan")
            in_isophote = np.zeros(ix.size, dtype=bool)
            in_isophote[on_image] = observed_region[iy[on_image], ix[on_image]]

            radius_sphere = SPHERICAL_RADIUS_KPC[band]
            method_masks = {
                "isophotal_particle_selection": in_isophote,
                "b010_linked_group": data["link_masks"][0.10],
                "b012_linked_group": data["link_masks"][0.12],
                "b020_linked_group": data["link_masks"][0.20],
                "projected_aperture_30_pkpc": radius2_pkpc <= 30.0,
                "Darwin_calibrated_sphere": data["radius3d_pkpc"] <= radius_sphere,
                "central_bound_stars_reference": data["central_mask"],
            }
            reference_mass = float(data["mass_msun"][data["central_mask"]].sum())
            reference_r50 = weighted_r50(radius2_pkpc, flux[:, iband], data["central_mask"])
            isophotal_mass = float(data["mass_msun"][in_isophote].sum())
            isophotal_r50 = weighted_r50(radius2_pkpc, flux[:, iband], in_isophote)
            for method, selected in method_masks.items():
                selected_mass = float(data["mass_msun"][selected].sum())
                selected_r50 = weighted_r50(radius2_pkpc, flux[:, iband], selected)
                selected_r50_mass = weighted_r50(radius2_pkpc, data["mass_msun"], selected)
                iou, completeness, precision = rv.light_scores(flux[:, iband], selected, in_isophote)
                rows.append({
                    **common,
                    "band": band,
                    "method": method,
                    "calibrated_sphere_radius_pkpc": radius_sphere if method == "Darwin_calibrated_sphere" else "",
                    "observed_isophote_area_arcsec2": observed_area,
                    "observed_isophote_pixel_count": int(observed_region.sum()),
                    "isophote_status": isophote_status,
                    "central_pixel_detected": int(central_pixel_detected),
                    "nearest_bright_pixel_fallback_used": int(
                        isophote_status == "nearest_bright_pixel_fallback"),
                    "isophote_touches_image_edge": int(edge_contact),
                    "isophote_centroid_offset_pixels": region_centroid_offset_pix,
                    "nearest_isophote_pixel_offset_pixels": nearest_region_pixel_offset_pix,
                    "selected_particle_count": int(selected.sum()),
                    "light_weighted_IoU": float(iou),
                    "light_completeness": float(completeness),
                    "light_precision": float(precision),
                    "stellar_mass_msun": selected_mass,
                    "central_bound_stellar_mass_msun": reference_mass,
                    "log10_mass_ratio_vs_central_bound_dex": (
                        float(np.log10(selected_mass / reference_mass))
                        if selected_mass > 0 and reference_mass > 0 else float("nan")),
                    "isophotal_stellar_mass_msun": isophotal_mass,
                    "log10_mass_ratio_vs_isophote_dex": (
                        float(np.log10(selected_mass / isophotal_mass))
                        if selected_mass > 0 and isophotal_mass > 0 else float("nan")),
                    "projected_R50_pkpc": selected_r50,
                    "projected_R50_mass_pkpc": selected_r50_mass,
                    "central_bound_R50_pkpc": reference_r50,
                    "central_bound_R50_mass_pkpc": weighted_r50(
                        radius2_pkpc, data["mass_msun"], data["central_mask"]),
                    "isophotal_R50_pkpc": isophotal_r50,
                    "log10_R50_ratio_vs_central_bound_dex": (
                        float(np.log10(selected_r50 / reference_r50))
                        if selected_r50 > 0 and reference_r50 > 0 else float("nan")),
                    "log10_R50_ratio_vs_isophote_dex": (
                        float(np.log10(selected_r50 / isophotal_r50))
                        if selected_r50 > 0 and isophotal_r50 > 0 else float("nan")),
                    "psf_fwhm_arcsec": float(sb.FWHM_AS),
                    "observer_redshift_is_counterfactual": int(abs(observer_redshift - z_phys) > 1.0e-5),
                })
    return rows


def run_rank(rank: int, output_dir: Path | None = None):
    output_dir = OUT if output_dir is None else output_dir
    output_dir.mkdir(parents=True, exist_ok=True)
    outpath = output_dir / f"tree_rank_{rank:02d}_growth_scores.csv"
    if outpath.exists():
        raise FileExistsError(f"Refusing to overwrite existing growth scores: {outpath}")
    rows = []
    grid_cache = {}
    for snapshot in TARGET_SNAPS:
        record = epoch_rows(rank)[snapshot]
        state = load_epoch(rank, snapshot, record)
        print(f"START rank={rank} snap={snapshot} z={state['z_phys']:.6f} "
              f"nstar={state['n_particles']} central={state['central_particle_count']} "
              f"catalog_mass_ratio={state['catalog_mass_ratio']:.6f}", flush=True)
        for z_observer in OBSERVER_REDSHIFTS:
            rows.extend(render_epoch(state, z_observer, grid_cache))
            print(f"  snap={snapshot} observer_z={z_observer:.6f} done", flush=True)
        del state
    if not rows:
        raise RuntimeError(f"No scores generated for rank {rank}")
    with outpath.open("x", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    print(f"DONE rank={rank} rows={len(rows)} output={outpath}", flush=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample-rank", type=int, required=True)
    parser.add_argument("--output-directory", type=Path, default=OUT)
    args = parser.parse_args()
    if args.sample_rank not in locked_cohort():
        raise SystemExit("sample rank is not in the 36-tree common-central cohort")
    run_rank(args.sample_rank, args.output_directory.resolve())


if __name__ == "__main__":
    main()
