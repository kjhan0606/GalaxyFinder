#!/usr/bin/env python3
"""Measure one CAMELS Swift-EAGLE halo with the locked PGalF image setup."""

from __future__ import annotations

import argparse
import csv
import json
import os
import sys
from pathlib import Path

os.environ.setdefault("SPS_HOME", "/home/kjhan/fsps_data")
os.environ.update(
    OMP_NUM_THREADS="1",
    OPENBLAS_NUM_THREADS="1",
    MKL_NUM_THREADS="1",
    NUMBA_NUM_THREADS="1",
)

import h5py
import numpy as np
from astropy.cosmology import FlatLambdaCDM
from scipy.signal import fftconvolve


PROJECT = Path(__file__).resolve().parents[1]
REVIEW = Path(__file__).resolve().parent
DATA_ROOT = Path("/gpfs/kjhan/CAMELS/Swift-EAGLE/L25n256/CV")
MANIFEST = REVIEW / "sample_manifest.csv"
RESULTS = REVIEW / "results"
FSPS_GRID = PROJECT / "PGalF_major_revision_20261013/inputs/supplement/fsps_z0.190.npz"
SNAPSHOT = 90
HUBBLE = 0.6711
OMEGA_M = 0.3
BOX_MPC = 37.25227241
OBSERVED_REDSHIFT = 0.19043796228481114
Z_SUN_FSPS = 0.0142
PIX_AS = 0.2
STAMP_AS = 180.0
FINE_B = np.array([0.02, 0.03, 0.05, 0.06, 0.07, 0.08, 0.09,
                   0.10, 0.11, 0.12, 0.13, 0.14, 0.15, 0.18, 0.20])
LOS = (("x", 1, 2), ("y", 2, 0), ("z", 0, 1))
BANDS = (("V", 0, 26.5), ("g", 1, 27.0))
FOF_UNASSIGNED = 2147483647
_SNAPSHOT_ID_CACHE: dict[int, tuple[np.ndarray, np.ndarray, np.ndarray]] = {}

sys.path.insert(0, str(PROJECT / "PGalF_opus_revision_20261013/scripts/source"))
import sb_vs_b as sb  # noqa: E402
import sb_revise as rv  # noqa: E402


def read_manifest_row(cv_box: int, group_index: int) -> dict[str, str]:
    with MANIFEST.open(newline="") as stream:
        for row in csv.DictReader(stream):
            if int(row["cv_box"]) == cv_box and int(row["fof_group_index"]) == group_index:
                return row
    raise KeyError(f"No selected halo CV_{cv_box}, group {group_index}")


def clear_snapshot_id_cache(cv_box: int | None = None) -> None:
    """Release cached full-snapshot membership arrays between CAMELS boxes."""
    if cv_box is None:
        _SNAPSHOT_ID_CACHE.clear()
    else:
        _SNAPSHOT_ID_CACHE.pop(cv_box, None)


def snapshot_id_index(cv_box: int, snapshot_path: Path):
    """Load full-snapshot IDs once per box for a batch of selected FoF groups."""
    if cv_box not in _SNAPSHOT_ID_CACHE:
        with h5py.File(snapshot_path, "r") as snap:
            stars = snap["PartType4"]
            particle_ids = stars["ParticleIDs"][:]
            fof_ids = stars["FOFGroupIDs"][:]
        sorted_particle_ids = np.sort(particle_ids)
        _SNAPSHOT_ID_CACHE[cv_box] = (particle_ids, fof_ids, sorted_particle_ids)
    return _SNAPSHOT_ID_CACHE[cv_box]


def read_input_stars(
    cv_box: int, group_index: int, expected_central: int, membership_source: str
):
    catalog_path = DATA_ROOT / f"CV_{cv_box}" / f"groups_{SNAPSHOT:03d}.hdf5"
    snapshot_path = Path("/gpfs/kjhan/CAMELS/Swift-EAGLE/L25n256/raw") / f"CV_{cv_box}" / f"snapshot_{SNAPSHOT:03d}.hdf5"

    with h5py.File(catalog_path, "r") as cat:
        group = cat["Group"]
        group_lengths = group["GroupLen"][:]
        group_star_catalog_count = int(group["GroupLenType"][group_index, 4])
        first_sub = int(group["GroupFirstSub"][group_index])
        n_sub = int(group["GroupNsubs"][group_index])
        sub_len_type = cat["Subhalo/SubhaloLenType"][:]
        sub_mass_type = cat["Subhalo/SubhaloMassType"][:]
        sub_cm = cat["Subhalo/SubhaloCM"][:]
        sub_pos = cat["Subhalo/SubhaloPos"][:]
        id_data = cat["IDs/ID"]
        group_id_start = int(np.sum(group_lengths[:group_index], dtype=np.int64))
        group_id_end = group_id_start + int(group_lengths[group_index])
        group_member_ids = np.asarray(id_data[group_id_start:group_id_end], dtype=np.int64)

        local_subs = np.arange(first_sub, first_sub + n_sub, dtype=np.int64)
        local_star_counts = sub_len_type[local_subs, 4]
        central_offset = int(np.argmax(local_star_counts))
        central_index = int(local_subs[central_offset])
        central_nstar = int(local_star_counts[central_offset])
        if central_index != expected_central:
            raise RuntimeError(
                f"Manifest/catalogue central mismatch: {expected_central} != {central_index}"
            )

        sub_id_offsets: dict[int, tuple[int, int]] = {}
        running = group_id_start
        for sub_index in local_subs:
            lengths = sub_len_type[int(sub_index)]
            n_this = int(np.sum(lengths, dtype=np.int64))
            sub_id_offsets[int(sub_index)] = (running, n_this)
            running += n_this

        if running > id_data.shape[0]:
            raise RuntimeError("Subhalo ID offsets exceed the catalogue ID array")
        if running > group_id_start + int(group_lengths[group_index]):
            raise RuntimeError("Subhalo ID offsets exceed the selected FOF group block")

        satellite_indices = [
            int(index) for index in local_subs
            if int(index) != central_index and int(sub_len_type[int(index), 4]) >= 100
        ]
        satellite_id_blocks = []
        for sub_index in satellite_indices:
            start, count = sub_id_offsets[sub_index]
            if count:
                satellite_id_blocks.append(np.asarray(id_data[start:start + count], dtype=np.int64))

        central_start, central_block_count = sub_id_offsets[central_index]
        central_id_block = np.asarray(
            id_data[central_start:central_start + central_block_count], dtype=np.int64
        )
        central_center_mpc = np.asarray(sub_cm[central_index], dtype=np.float64) / (1000.0 * HUBBLE)
        central_potential_mpc = np.asarray(sub_pos[central_index], dtype=np.float64) / (1000.0 * HUBBLE)
        potential_delta_mpc = central_potential_mpc - central_center_mpc
        potential_delta_mpc -= BOX_MPC * np.round(potential_delta_mpc / BOX_MPC)
        central_mass_catalog = float(sub_mass_type[central_index, 4])

    particle_ids, fof_ids, sorted_particle_ids = snapshot_id_index(cv_box, snapshot_path)
    with h5py.File(snapshot_path, "r") as snap:
        stars = snap["PartType4"]

        def star_ids_from_subhalo_block(block: np.ndarray, expected: int, label: str) -> np.ndarray:
            locations = np.searchsorted(sorted_particle_ids, block)
            found = locations < sorted_particle_ids.size
            valid_locations = np.flatnonzero(found)
            found[valid_locations] &= sorted_particle_ids[locations[valid_locations]] == block[valid_locations]
            result = block[found]
            if result.size != expected:
                raise RuntimeError(
                    f"Subhalo {label} ID block contains {result.size} snapshot stars; expected {expected}"
                )
            if np.unique(result).size != result.size:
                raise RuntimeError(f"Duplicate star IDs in subhalo {label} ID block")
            return result

        central_ids = star_ids_from_subhalo_block(
            central_id_block, central_nstar, str(central_index)
        )
        satellite_id_chunks = [
            star_ids_from_subhalo_block(
                block, int(sub_len_type[si, 4]), str(si)
            )
            for si, block in zip(satellite_indices, satellite_id_blocks)
        ]
        satellite_ids = (
            np.unique(np.concatenate(satellite_id_chunks))
            if satellite_id_chunks else np.empty(0, dtype=np.int64)
        )
        fof_mask = fof_ids == group_index + 1
        if not np.any(fof_mask):
            raise RuntimeError(f"No raw star particles have FOFGroupID={group_index + 1}")

        raw_group_ids = particle_ids[fof_mask]
        catalogue_member_mask = np.isin(particle_ids, group_member_ids, assume_unique=False)
        catalogue_group_star_ids = particle_ids[catalogue_member_mask]
        if catalogue_group_star_ids.size != group_star_catalog_count:
            raise RuntimeError(
                "Catalogue group member IDs do not reproduce GroupLenType star count: "
                f"{catalogue_group_star_ids.size} != {group_star_catalog_count}"
            )
        if membership_source == "snapshot_fof":
            selected_member_mask = fof_mask
        elif membership_source == "catalogue_ids":
            selected_member_mask = catalogue_member_mask
        else:
            raise ValueError(f"Unknown membership source: {membership_source}")
        selected_member_ids = particle_ids[selected_member_mask]
        satellite_in_selected = np.isin(
            selected_member_ids, satellite_ids, assume_unique=False
        )
        input_mask = selected_member_mask.copy()
        input_mask[np.flatnonzero(selected_member_mask)] &= ~satellite_in_selected

        input_indices = np.flatnonzero(input_mask)
        coords_mpc = stars["Coordinates"][input_indices].astype(np.float64)
        current_mass_msun = stars["Masses"][input_indices].astype(np.float64) * 1.0e10
        formed_mass_msun = stars["InitialMasses"][input_indices].astype(np.float64) * 1.0e10
        birth_scale = stars["BirthScaleFactors"][input_indices].astype(np.float64)
        metallicity = stars["MetalMassFractions"][input_indices].astype(np.float64)
        selected_ids = particle_ids[input_indices]
        input_fof_ids = fof_ids[input_indices]

    if np.unique(selected_ids).size != selected_ids.size:
        raise RuntimeError("Duplicate stellar IDs in the selected input list")
    if membership_source == "snapshot_fof" and not np.all(input_fof_ids == group_index + 1):
        raise RuntimeError("A selected star has an unexpected raw FOFGroupID")

    # Center the image, fixed aperture, and radial measurements on the
    # central subhalo's minimum-potential position. FoF membership is
    # translation invariant, so this changes only the spatial comparison.
    delta = coords_mpc - central_potential_mpc
    delta -= BOX_MPC * np.round(delta / BOX_MPC)
    coords_mpc = delta
    coords_hinv_mpc = delta * HUBBLE

    central_match = int(np.isin(central_ids, raw_group_ids).sum())
    central_present = np.isin(central_ids, particle_ids)
    central_fof_ids = fof_ids[np.isin(particle_ids, central_ids)]
    central_unassigned = int(np.count_nonzero(central_fof_ids == FOF_UNASSIGNED))
    central_other_fof = int(np.count_nonzero(
        (central_fof_ids != group_index + 1)
        & (central_fof_ids != FOF_UNASSIGNED)
    ))
    satellite_match = int(np.isin(raw_group_ids, satellite_ids, assume_unique=False).sum())
    metadata = {
        "membership_source": membership_source,
        "catalogue_central_star_count": central_nstar,
        "catalogue_central_stellar_mass_1e10_msun_per_h": central_mass_catalog,
        "raw_fof_star_count": int(fof_mask.sum()),
        "catalogue_fof_star_count": group_star_catalog_count,
        "catalogue_group_star_ids_matching_snapshot": int(catalogue_group_star_ids.size),
        "catalogue_group_star_ids_matching_raw_fof": int(
            np.count_nonzero(np.isin(catalogue_group_star_ids, raw_group_ids))
        ),
        "raw_fof_star_ids_matching_catalogue_group": int(
            np.count_nonzero(np.isin(raw_group_ids, catalogue_group_star_ids))
        ),
        "selected_membership_star_count_before_satellite_mask": int(selected_member_ids.size),
        "resolved_satellite_ids_in_selected_membership": int(satellite_in_selected.sum()),
        "central_catalogue_ids_matching_raw_fof": central_match,
        "central_catalogue_id_count": int(central_ids.size),
        "central_catalogue_id_match_fraction": central_match / central_ids.size if central_ids.size else 0.0,
        "central_catalogue_ids_found_in_snapshot": int(central_present.sum()),
        "central_catalogue_ids_unassigned_to_fof": central_unassigned,
        "central_catalogue_ids_assigned_to_other_fof": central_other_fof,
        "subhalocm_to_potential_offset_kpc": float(np.linalg.norm(potential_delta_mpc) * 1000.0),
        "potential_delta_from_subhalocm_mpc": potential_delta_mpc.tolist(),
        "radial_measurement_centre": "central subhalo minimum-potential position",
        "resolved_satellite_subhalo_count": len(satellite_indices),
        "resolved_satellite_catalogue_star_count": int(satellite_ids.size),
        "resolved_satellite_ids_matching_raw_fof": satellite_match,
        "input_star_count_after_satellite_mask": int(input_indices.size),
        "mass_units_source": "PartType4 Masses / InitialMasses: 1e10 Msun per unit; h exponent 0",
        "position_units_source": "PartType4 Coordinates: comoving Mpc; h exponent 0",
    }
    return coords_mpc, coords_hinv_mpc, current_mass_msun * HUBBLE, formed_mass_msun, birth_scale, metallicity, metadata


def star_fluxes(formed_mass_msun, birth_scale, metallicity, grids, camels_cosmo, fsps_cosmo):
    age_grid, logz_grid, unit_flux = grids
    valid_a = np.clip(birth_scale, 1.0e-4, 1.0)
    birth_redshift = 1.0 / valid_a - 1.0
    age = camels_cosmo.age(0.0).value - camels_cosmo.age(birth_redshift).value
    age = np.clip(age, age_grid[0], age_grid[-1])
    logz = np.log10(np.clip(metallicity, 1.0e-5, 0.05) / Z_SUN_FSPS)

    log_age_grid = np.log(age_grid)
    x = np.interp(np.log(age), log_age_grid, np.arange(age_grid.size))
    y = np.interp(logz, logz_grid, np.arange(logz_grid.size))
    x0 = np.clip(np.floor(x).astype(np.int64), 0, age_grid.size - 2)
    y0 = np.clip(np.floor(y).astype(np.int64), 0, logz_grid.size - 2)
    tx = x - x0
    ty = y - y0

    flux = np.empty((age.size, 2), dtype=np.float64)
    distance_factor = (
        fsps_cosmo.luminosity_distance(OBSERVED_REDSHIFT).value
        / camels_cosmo.luminosity_distance(OBSERVED_REDSHIFT).value
    ) ** 2
    for band in range(2):
        table = np.log(unit_flux[:, :, band])
        f00 = table[x0, y0]
        f10 = table[x0 + 1, y0]
        f01 = table[x0, y0 + 1]
        f11 = table[x0 + 1, y0 + 1]
        flux[:, band] = formed_mass_msun * np.exp(
            (1 - tx) * (1 - ty) * f00
            + tx * (1 - ty) * f10
            + (1 - tx) * ty * f01
            + tx * ty * f11
        ) * distance_factor
    return flux, age, distance_factor


def weighted_r50(radius_kpc, weights, selection):
    use = selection & np.isfinite(weights) & (weights > 0)
    if not np.any(use):
        return float("nan")
    radius = radius_kpc[use]
    weight = weights[use]
    order = np.argsort(radius)
    csum = np.cumsum(weight[order])
    return float(radius[order[np.searchsorted(csum, 0.5 * csum[-1], side="left")]])


def score(selection, target, weight, radius_kpc):
    lg = float(weight[selection].sum())
    lt = float(weight[target].sum())
    both = selection & target
    lb = float(weight[both].sum())
    union = lg + lt - lb
    return {
        "iou": lb / union if union > 0 else float("nan"),
        "completeness": lb / lt if lt > 0 else float("nan"),
        "precision": lb / lg if lg > 0 else float("nan"),
        "L_selected": lg,
        "L_contour": lt,
        "L_common": lb,
        "luminosity_ratio_minus_one": lg / lt - 1 if lt > 0 else float("nan"),
        "R50_selected_kpc": weighted_r50(radius_kpc, weight, selection),
        "R50_contour_kpc": weighted_r50(radius_kpc, weight, target),
    }


def measure(
    cv_box: int,
    group_index: int,
    save_images: bool,
    tag: str = "",
    membership_source: str = "catalogue_ids",
) -> Path:
    suffix = f"_{tag}" if tag else ""
    result_dir = RESULTS / f"CV{cv_box:02d}_group{group_index:05d}{suffix}"
    if result_dir.exists():
        raise FileExistsError(f"Refusing to overwrite {result_dir}")
    manifest = read_manifest_row(cv_box, group_index)
    expected_central = int(manifest["central_subhalo_index"])
    coords_mpc, coords_hinv_mpc, mass_hinv, formed_mass_msun, birth_scale, metallicity, data_meta = read_input_stars(
        cv_box, group_index, expected_central, membership_source
    )

    camels_cosmo = FlatLambdaCDM(H0=100.0 * HUBBLE, Om0=OMEGA_M)
    fsps_cosmo = FlatLambdaCDM(H0=72.0, Om0=0.27)
    with np.load(FSPS_GRID) as grid:
        grid_redshift = float(grid["z"])
        age_grid = grid["age"]
        logz_grid = grid["logz"]
        unit_flux = grid["flux"]
    if not np.isclose(grid_redshift, OBSERVED_REDSHIFT, rtol=0.0, atol=1e-12):
        raise RuntimeError(f"Unexpected FSPS grid redshift: {grid_redshift}")
    flux, ages_gyr, distance_factor = star_fluxes(
        formed_mass_msun, birth_scale, metallicity,
        (age_grid, logz_grid, unit_flux), camels_cosmo, fsps_cosmo,
    )

    da_mpc = camels_cosmo.angular_diameter_distance(OBSERVED_REDSHIFT).value
    arcsec_per_mpc = 206265.0 / da_mpc
    n = int(round(2 * STAMP_AS / PIX_AS))
    pixels_per_mpc = arcsec_per_mpc / PIX_AS
    projected = []
    for los_name, axis_x, axis_y in LOS:
        ix = np.rint(coords_mpc[:, axis_x] * pixels_per_mpc).astype(np.int32) + n // 2
        iy = np.rint(coords_mpc[:, axis_y] * pixels_per_mpc).astype(np.int32) + n // 2
        on = (ix >= 0) & (ix < n) & (iy >= 0) & (iy < n)
        radius_kpc = np.hypot(coords_mpc[:, axis_x], coords_mpc[:, axis_y]) * 1000.0
        aperture = radius_kpc <= 30.0
        projected.append((los_name, axis_x, axis_y, ix, iy, on, radius_kpc, aperture))

    seed = sb.nucleus_index(*[coords_hinv_mpc[:, axis] for axis in range(3)])
    cm_delta_from_potential = -np.asarray(
        data_meta["potential_delta_from_subhalocm_mpc"], dtype=np.float64
    )
    seed_to_subhalocm = coords_mpc[seed] - cm_delta_from_potential
    seed_to_subhalocm -= BOX_MPC * np.round(seed_to_subhalocm / BOX_MPC)
    rows: list[dict[str, object]] = []
    image_store: dict[str, np.ndarray] = {}
    psf_kernels = [
        sb.psf_kernel(PIX_AS, sb.FWHM_AS, sb.MOFFAT_BETA, sb.WING_A[index], sb.KERNEL_AS)[0]
        for index in range(2)
    ]

    targets = {}
    for view_index, item in enumerate(projected):
        los_name, axis_x, axis_y, ix, iy, on, radius_kpc, aperture = item
        seed_offset_arcsec = float(
            np.hypot(coords_mpc[seed, axis_x], coords_mpc[seed, axis_y]) * arcsec_per_mpc
        )
        for band_name, band_index, threshold in BANDS:
            image = sb.bin_flux(ix, iy, flux[:, band_index], n)
            image = fftconvolve(image, psf_kernels[band_index], mode="same")
            image[image < 0] = 0.0
            mu, mask, labels, nlabel = rv.mu_labels(image, threshold)
            peak_y, peak_x = np.unravel_index(int(np.argmax(image)), image.shape)
            peak_offset_arcsec = float(
                np.hypot(peak_x - n // 2, peak_y - n // 2) * PIX_AS
            )
            if save_images:
                image_store[f"image_{los_name}_{band_name}"] = image.astype(np.float32)
            for variant, py, px in (
                (
                    "subhalocm",
                    int(np.rint(cm_delta_from_potential[axis_y] * pixels_per_mpc)) + n // 2,
                    int(np.rint(cm_delta_from_potential[axis_x] * pixels_per_mpc)) + n // 2,
                ),
                ("potential", n // 2, n // 2),
                ("light_peak", int(peak_y), int(peak_x)),
            ):
                component, fallback, _clipped = rv.component_at(labels, mask, py, px)
                region, area = rv.region_of(labels, component, PIX_AS)
                target = rv.sample_inside(region, ix, iy, on)
                targets[(view_index, band_index, variant)] = target
                if save_images:
                    image_store[f"mask_{variant}_{los_name}_{band_name}"] = region.astype(np.uint8)

                rows.append({
                    "cv_box": cv_box,
                    "fof_group_index": group_index,
                    "line_of_sight": los_name,
                    "band": band_name,
                    "target_variant": variant,
                    "method": "aperture30pkpc",
                    "b_star": "",
                    **score(aperture, target, flux[:, band_index], radius_kpc),
                    "target_area_arcsec2": area,
                    "target_npixel": int(region.sum()),
                    "target_nstar": int(target.sum()),
                    "target_component_fallback": int(fallback),
                    "density_seed_offset_arcsec": seed_offset_arcsec,
                    "light_peak_offset_arcsec": peak_offset_arcsec,
                    "n_components": nlabel,
                })
            del image, mu, mask, labels, region

    for b_star in FINE_B:
        link = b_star * (mass_hinv / (2.7755e11 * OMEGA_M)) ** (1.0 / 3.0)
        roots = sb.fof_largest(
            coords_hinv_mpc[:, 0], coords_hinv_mpc[:, 1], coords_hinv_mpc[:, 2], link
        )
        selected = sb.peak_mask(roots, *[coords_hinv_mpc[:, axis] for axis in range(3)], seed)
        for view_index, item in enumerate(projected):
            los_name, _, _, _, _, _, radius_kpc, _ = item
            seed_offset_arcsec = float(
                np.hypot(coords_mpc[seed, item[1]], coords_mpc[seed, item[2]]) * arcsec_per_mpc
            )
            for band_name, band_index, _ in BANDS:
                for variant in ("subhalocm", "potential", "light_peak"):
                    target = targets[(view_index, band_index, variant)]
                    peak_offset_arcsec = next(
                        float(r["light_peak_offset_arcsec"])
                        for r in rows
                        if r["line_of_sight"] == los_name
                        and r["band"] == band_name
                        and r["target_variant"] == variant
                    )
                    rows.append({
                        "cv_box": cv_box,
                        "fof_group_index": group_index,
                        "line_of_sight": los_name,
                        "band": band_name,
                        "target_variant": variant,
                        "method": "FoF",
                        "b_star": float(b_star),
                        **score(selected, target, flux[:, band_index], radius_kpc),
                        "target_area_arcsec2": "",
                        "target_npixel": "",
                        "target_nstar": int(target.sum()),
                        "target_component_fallback": "",
                        "density_seed_offset_arcsec": seed_offset_arcsec,
                        "light_peak_offset_arcsec": peak_offset_arcsec,
                        "n_components": "",
                    })

    result_dir.mkdir(parents=True)
    metrics_path = result_dir / "metrics.csv"
    with metrics_path.open("x", newline="") as stream:
        fields = list(rows[0].keys())
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)

    meta = {
        "cv_box": cv_box,
        "snapshot": SNAPSHOT,
        "simulation_redshift": float(manifest["redshift"]),
        "mock_observed_redshift": OBSERVED_REDSHIFT,
        "mock_redshift_note": "The z=0 structure and stellar ages are placed at z_obs=0.1904379623 for matched image scale and redshifted FSPS photometry; this is not a z=0.19 simulation output.",
        "fof_group_index": group_index,
        "fof_group_id_raw": int(manifest["fof_group_id_raw"]),
        "fof_mass_1e10_msun_per_h": float(manifest["fof_mass_1e10_msun_per_h"]),
        "central_subhalo_index": expected_central,
        "input_count": int(coords_mpc.shape[0]),
        "density_seed_index": int(seed),
        "density_seed_offset_from_subhalo_cm_kpc": float(np.linalg.norm(seed_to_subhalocm) * 1000.0),
        "density_seed_offset_from_subhalo_potential_kpc": float(np.linalg.norm(coords_mpc[seed]) * 1000.0),
        "target_variants": ["subhalocm", "potential", "light_peak"],
        "minimum_age_gyr": float(np.min(ages_gyr)),
        "median_age_gyr": float(np.median(ages_gyr)),
        "maximum_age_gyr": float(np.max(ages_gyr)),
        "median_metal_mass_fraction": float(np.median(metallicity)),
        "fsps_grid_read_only_path": str(FSPS_GRID),
        "fsps_version": "0.4.7",
        "distance_flux_factor_fsps_over_camels_squared": float(distance_factor),
        "angular_diameter_distance_mpc": float(da_mpc),
        "pixel_scale_arcsec": PIX_AS,
        "stamp_half_width_arcsec": STAMP_AS,
        "b_grid": FINE_B.tolist(),
        **data_meta,
    }
    (result_dir / "metadata.json").write_text(json.dumps(meta, indent=2) + "\n")
    if save_images:
        np.savez_compressed(result_dir / "images_and_contours.npz", **image_store)
    print(f"Wrote {len(rows)} metric rows to {metrics_path}")
    print(json.dumps({k: meta[k] for k in (
        "input_count", "raw_fof_star_count", "catalogue_central_star_count",
        "resolved_satellite_subhalo_count", "resolved_satellite_ids_matching_raw_fof",
        "distance_flux_factor_fsps_over_camels_squared"
    )}, indent=2))
    return result_dir


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cv-box", type=int, required=True)
    parser.add_argument("--group-index", type=int, required=True)
    parser.add_argument("--save-images", action="store_true")
    parser.add_argument("--tag", default="")
    parser.add_argument(
        "--membership-source",
        choices=("snapshot_fof", "catalogue_ids"),
        default="catalogue_ids",
    )
    args = parser.parse_args()
    measure(
        args.cv_box, args.group_index, args.save_images, args.tag,
        args.membership_source,
    )


if __name__ == "__main__":
    main()
