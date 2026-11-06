#!/usr/bin/env python3
"""Track central-star assembly and stellar-population size weighting in TNG300-1.

The script uses the locked 36-branch TNG cohort and the exact particle
membership cutouts used by the current growth analysis. It writes a closed
particle-mass ledger, selection capture fractions by assembly class, radial
stellar-population profiles, and controlled half-light-radius growth
measurements. Branches, rather than projections or particle records, are the
statistical units.
"""
from __future__ import annotations

import csv
import hashlib
import json
import math
import sys
from collections import defaultdict
from pathlib import Path

import h5py
import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import fftconvolve

BASE = Path(__file__).resolve().parents[1]
OUT = Path(__file__).resolve().parent
OUTPUT_DIR = OUT / "results"
sys.path.insert(0, str(BASE))
import render_tng_growth_crossrender as renderer  # noqa: E402

sb = renderer.sb
H = renderer.H_TNG
LOS_AXES = renderer.AXES
BANDS = ("V", "g")
SNAPS = (50, 67, 84)
ENDPOINT_INTERVALS = ((50, 67), (67, 84), (50, 84))
OBSERVER_FIXED = renderer.OBSERVER_REDSHIFTS[0]
METHODS = (
    "isophotal_particle_selection", "b010_linked_group", "b012_linked_group",
    "b020_linked_group", "projected_aperture_30_pkpc", "Darwin_calibrated_sphere",
)
BOOTSTRAPS = 20000
SEED = 20261106
RADIAL_EDGES = (0.0, 0.5, 1.0, 2.0, 4.0, 8.0, np.inf)
RADIAL_LABELS = ("0-0.5", "0.5-1", "1-2", "2-4", "4-8", ">8")


def read_csv(path: Path):
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream))


def write_csv(path: Path, rows):
    if not rows:
        raise RuntimeError(f"No rows to write: {path}")
    with path.open("x", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def write_json(path: Path, value):
    with path.open("x") as stream:
        json.dump(value, stream, indent=2, sort_keys=True)
        stream.write("\n")


def epoch_index():
    result = defaultdict(dict)
    for row in read_csv(BASE / "lineage_epochs.csv"):
        rank, snap = int(row["sample_rank"]), int(row["snapshot"])
        if snap in SNAPS:
            result[rank][snap] = row
    return result


def cohort_ranks():
    return sorted(int(row["sample_rank"]) for row in read_csv(BASE / "lineage_cohort.csv")
                  if int(row["central_at_all_targets"]) == 1)


def membership_path(rank: int, snap: int, row: dict[str, str]) -> Path:
    group = int(row["progenitor_fof_group_id"])
    return BASE / "tree_epoch_membership" / f"snap{snap:03d}_rank{rank:02d}_group{group:06d}.h5"


def load_epoch(rank: int, snap: int, row: dict[str, str]):
    """Load both the full parent-FoF particle table and renderer-matched stars."""
    path = membership_path(rank, snap, row)
    with h5py.File(path, "r") as handle:
        p = handle["PartType4"]
        coords = p["Coordinates"][:].astype(np.float64)
        masses = p["Masses"][:].astype(np.float64)
        initial = p["GFM_InitialMass"][:].astype(np.float64)
        metal = p["GFM_Metallicity"][:].astype(np.float64)
        birth = p["GFM_StellarFormationTime"][:].astype(np.float64)
        ids = p["ParticleIDs"][:].astype(np.uint64)
        central = p["CentralSubhaloParticle"][:].astype(bool)
        satellite = p["RemoveSatelliteGalaxy"][:].astype(bool)
        z = float(handle.attrs["Redshift"])
        box_ckpc_h = float(handle.attrs.get("BoxSize", renderer.BOX_CKPC_H))
        group_id = int(handle.attrs["GroupID"])
        central_id = int(handle.attrs["CentralSubhaloID"])
    if group_id != int(row["progenitor_fof_group_id"]):
        raise RuntimeError(f"FoF group mismatch rank={rank} snapshot={snap}")
    if central_id != int(row["progenitor_subhalo_id"]):
        raise RuntimeError(f"SubLink central mismatch rank={rank} snapshot={snap}")
    if len(np.unique(ids)) != ids.size:
        raise RuntimeError(f"Duplicate ParticleIDs within {path.name}")
    if np.any(central & satellite):
        raise RuntimeError(f"Central and satellite masks overlap in {path.name}")

    a_snap = 1.0 / (1.0 + z)
    valid = np.isfinite(birth) & (birth > 0.0) & (birth <= a_snap)
    central_valid = central & valid
    expected_central_n = int(row["stellar_particle_count"])
    if int(central.sum()) != expected_central_n:
        raise RuntimeError(
            f"Central count mismatch rank={rank} snap={snap}: "
            f"{central.sum()} != {expected_central_n}")
    if np.any(central_valid & satellite):
        raise RuntimeError(f"A central member is masked as a satellite in {path.name}")

    center = np.array([float(row[f"subhalo_{axis}_ckpc_h"]) for axis in "xyz"], dtype=float)
    use_idx = np.flatnonzero(valid & ~satellite)
    delta = coords[use_idx] / 1000.0 - center / 1000.0
    box = box_ckpc_h / 1000.0
    delta -= box * np.round(delta / box)
    mass_h = masses[use_idx] * 1.0e10
    mass_msun = mass_h / H
    formed_mass_msun = initial[use_idx] * 1.0e10 / H
    birth_used = np.clip(birth[use_idx], 0.02, a_snap)
    a_grid, tau_grid, t_cum = renderer.time_grids(a_snap)
    birth_tau = np.interp(birth_used, a_grid, tau_grid)
    age = sb.proper_age(birth_used, a_snap, a_grid, t_cum)
    r3d_pkpc = np.linalg.norm(delta, axis=1) * a_snap / H * 1000.0

    seed = sb.nucleus_index(delta[:, 0], delta[:, 1], delta[:, 2])
    link_masks = {}
    for b in (0.10, 0.12, 0.20):
        link = b * np.cbrt(mass_h / (2.7755e11 * renderer.OM_TNG))
        labels = sb.fof_largest(delta[:, 0], delta[:, 1], delta[:, 2], link)
        link_masks[b] = labels == labels[seed]

    raw = {
        "path": path, "ids": ids, "coords": coords, "mass_code": masses,
        "initial_code": initial, "metal": metal, "birth_a": birth,
        "central": central, "satellite": satellite, "valid": valid,
        "central_valid": central_valid, "z": z, "a": a_snap,
        "box_ckpc_h": box_ckpc_h,
    }
    data = {
        "rank": rank, "snapshot": snap, "row": row, "z_phys": z,
        "a_snap": a_snap, "ids": ids[use_idx], "use_idx": use_idx,
        "delta": delta, "radius3d_pkpc": r3d_pkpc,
        "mass_msun": mass_msun, "mass_h": mass_h,
        "initial_msun": formed_mass_msun, "metal": metal[use_idx],
        "birth_a": birth_used, "birth_tau": birth_tau, "age_gyr": age,
        "central": central[use_idx], "link_masks": link_masks,
        "a_grid": a_grid, "tau_grid": tau_grid, "t_cum": t_cum,
        "n_input_stars": int(use_idx.size), "central_n": expected_central_n,
        "box_ckpc_h": box_ckpc_h,
    }
    cat_mass = float(row["stellar_mass_1e10_msun_h"]) * 1.0e10 / H
    particle_mass = float(mass_msun[data["central"]].sum())
    data["catalog_mass_ratio"] = particle_mass / cat_mass
    if abs(data["catalog_mass_ratio"] - 1.0) > 2e-4:
        raise RuntimeError(f"Central particle/catalogue mass mismatch rank={rank} snap={snap}")
    return raw, data


def id_lookup(reference_ids: np.ndarray, query_ids: np.ndarray):
    order = np.argsort(reference_ids)
    sorted_ids = reference_ids[order]
    pos = np.searchsorted(sorted_ids, query_ids)
    in_range = pos < sorted_ids.size
    found = np.zeros(query_ids.size, dtype=bool)
    found[in_range] = sorted_ids[pos[in_range]] == query_ids[in_range]
    ref_idx = np.full(query_ids.size, -1, dtype=np.int64)
    ref_idx[found] = order[pos[found]]
    return ref_idx, found


def assembly_transition(start: dict, end: dict, start_snap: int, end_snap: int, rank: int):
    start_c = np.flatnonzero(start["central_valid"])
    end_c = np.flatnonzero(end["central_valid"])
    start_ids = start["ids"][start_c]
    end_ids = end["ids"][end_c]
    end_to_start_c, end_found_start_c = id_lookup(start_ids, end_ids)
    persistent_end = end_found_start_c
    acquired_end = ~persistent_end
    born_after_start = end["birth_a"][end_c] > start["a"] + 1.0e-7
    new_end = acquired_end & born_after_start
    old_acquired_end = acquired_end & ~born_after_start

    start_to_end_c, start_found_end_c = id_lookup(end_ids, start_ids)
    lost_start = ~start_found_end_c

    start_mass = start["mass_code"][start_c] * 1.0e10 / H
    end_mass = end["mass_code"][end_c] * 1.0e10 / H
    persistent_delta = float(np.sum(
        end_mass[persistent_end] - start_mass[end_to_start_c[persistent_end]]))
    new_mass = float(end_mass[new_end].sum())
    old_acquired_mass = float(end_mass[old_acquired_end].sum())
    lost_mass = float(start_mass[lost_start].sum())
    start_total = float(start_mass.sum())
    end_total = float(end_mass.sum())
    delta = end_total - start_total
    reconstructed = persistent_delta + new_mass + old_acquired_mass - lost_mass

    # Source of pre-existing stars entering the central subhalo, from the
    # earlier parent-FoF cutout and its saved resolved-satellite mask.
    acquired_ids = end_ids[old_acquired_end]
    acquired_idx_end = end_c[old_acquired_end]
    start_idx_for_acquired, acquired_present = id_lookup(start["ids"], acquired_ids)
    acquired_sat = np.zeros(acquired_ids.size, dtype=bool)
    acquired_other = np.zeros(acquired_ids.size, dtype=bool)
    acquired_absent = ~acquired_present
    acquired_sat[acquired_present] = start["satellite"][start_idx_for_acquired[acquired_present]]
    acquired_other[acquired_present] = (
        ~start["satellite"][start_idx_for_acquired[acquired_present]]
        & ~start["central"][start_idx_for_acquired[acquired_present]])
    # A particle labelled as an earlier central member but absent from the
    # valid central set is retained as a distinct audit category.
    acquired_start_c_invalid = np.zeros(acquired_ids.size, dtype=bool)
    acquired_start_c_invalid[acquired_present] = (
        start["central"][start_idx_for_acquired[acquired_present]]
        & ~start["valid"][start_idx_for_acquired[acquired_present]])

    lost_ids = start_ids[lost_start]
    lost_idx_start = start_c[lost_start]
    end_idx_for_lost, lost_present_end = id_lookup(end["ids"], lost_ids)
    lost_sat_end = np.zeros(lost_ids.size, dtype=bool)
    lost_other_end = np.zeros(lost_ids.size, dtype=bool)
    lost_absent_end = ~lost_present_end
    lost_sat_end[lost_present_end] = end["satellite"][end_idx_for_lost[lost_present_end]]
    lost_other_end[lost_present_end] = (
        ~end["satellite"][end_idx_for_lost[lost_present_end]]
        & ~end["central"][end_idx_for_lost[lost_present_end]])

    result = {
        "sample_rank": rank, "snapshot_start": start_snap, "snapshot_end": end_snap,
        "redshift_start": start["z"], "redshift_end": end["z"],
        "n_start_central": int(start_c.size), "n_end_central": int(end_c.size),
        "n_persistent_central": int(persistent_end.sum()),
        "n_acquired_old": int(old_acquired_end.sum()),
        "n_newly_formed": int(new_end.sum()), "n_lost_central": int(lost_start.sum()),
        "mass_start_msun": start_total, "mass_end_msun": end_total,
        "delta_central_mass_msun": delta,
        "persistent_particle_mass_change_msun": persistent_delta,
        "acquired_preexisting_mass_msun": old_acquired_mass,
        "newly_formed_particle_mass_msun": new_mass,
        "lost_central_particle_mass_msun": lost_mass,
        "reconstructed_delta_msun": reconstructed,
        "closure_residual_msun": delta - reconstructed,
        "closure_residual_fraction_of_start": (delta - reconstructed) / start_total,
        "delta_over_start": delta / start_total,
        "persistent_change_over_start": persistent_delta / start_total,
        "acquired_old_over_start": old_acquired_mass / start_total,
        "new_formed_over_start": new_mass / start_total,
        "lost_over_start": lost_mass / start_total,
        "n_acquired_old_in_start_satellite_mask": int(acquired_sat.sum()),
        "n_acquired_old_in_start_parent_fof_other": int(acquired_other.sum()),
        "n_acquired_old_absent_from_start_parent_fof": int(acquired_absent.sum()),
        "n_acquired_old_start_central_invalid_birth": int(acquired_start_c_invalid.sum()),
        "mass_acquired_old_start_satellite_mask_msun": float(
            end["mass_code"][acquired_idx_end[acquired_sat]].sum() * 1.0e10 / H),
        "mass_acquired_old_start_parent_fof_other_msun": float(
            end["mass_code"][acquired_idx_end[acquired_other]].sum() * 1.0e10 / H),
        "mass_acquired_old_absent_start_parent_fof_msun": float(
            end["mass_code"][acquired_idx_end[acquired_absent]].sum() * 1.0e10 / H),
        "n_lost_present_in_end_satellite_mask": int(lost_sat_end.sum()),
        "n_lost_present_in_end_parent_fof_other": int(lost_other_end.sum()),
        "n_lost_absent_from_end_parent_fof": int(lost_absent_end.sum()),
    }
    class_labels = np.full(end["ids"].size, "not_central", dtype=object)
    # For end-central stars, distinguish persistence, newly formed stars,
    # and previously formed stars newly assigned to the central subhalo.
    central_classes = np.full(end_c.size, "unclassified_central", dtype=object)
    central_classes[persistent_end] = "persistent_since_start"
    central_classes[new_end] = "formed_since_start"
    old_positions = np.flatnonzero(old_acquired_end)
    central_classes[old_positions[acquired_sat]] = "acquired_old_from_resolved_satellite"
    central_classes[old_positions[acquired_other]] = "acquired_old_parent_fof_noncentral"
    central_classes[old_positions[acquired_absent]] = "acquired_old_absent_from_parent_fof"
    central_classes[old_positions[acquired_start_c_invalid]] = "acquired_old_invalid_start_membership"
    class_labels[end_c] = central_classes
    # Retain the exact source labels and end-central indices for population summaries.
    transition = {
        "class_labels_raw": class_labels,
        "end_central_indices": end_c,
        "persistent_end_ids": end_ids[persistent_end],
        "new_end_ids": end_ids[new_end],
        "old_acquired_end_ids": end_ids[old_acquired_end],
    }
    return result, transition


def mass_weighted_mean(values, weights):
    good = np.isfinite(values) & np.isfinite(weights) & (weights > 0)
    if not np.any(good):
        return float("nan")
    return float(np.average(values[good], weights=weights[good]))


def weighted_r50(radius, weights, selection):
    use = selection & np.isfinite(radius) & np.isfinite(weights) & (weights > 0.0)
    if not np.any(use):
        return float("nan")
    r, w = radius[use], weights[use]
    order = np.argsort(r)
    cumulative = np.cumsum(w[order])
    return float(r[order[np.searchsorted(cumulative, 0.5 * cumulative[-1], side="left")]])


def particle_flux(age, metal, formed_mass, observer_z, grid):
    """Bilinear interpolation matching sb.star_flux for supplied ages and Z."""
    age_grid, logz_grid, flux_grid = grid
    age = np.clip(np.asarray(age, dtype=float), age_grid[0], age_grid[-1])
    metal = np.clip(np.asarray(metal, dtype=float), 1.0e-5, 0.05)
    zsun = 0.0142
    lx, ly = np.log(age_grid), logz_grid
    x = np.interp(np.log(age), lx, np.arange(age_grid.size))
    y = np.interp(np.log10(metal / zsun), ly, np.arange(logz_grid.size))
    x0 = np.clip(np.floor(x).astype(int), 0, age_grid.size - 2)
    y0 = np.clip(np.floor(y).astype(int), 0, logz_grid.size - 2)
    tx, ty = x - x0, y - y0
    out = np.empty((age.size, flux_grid.shape[-1]), dtype=float)
    for band in range(flux_grid.shape[-1]):
        f = np.log(flux_grid[:, :, band])
        f00, f10 = f[x0, y0], f[x0 + 1, y0]
        f01, f11 = f[x0, y0 + 1], f[x0 + 1, y0 + 1]
        out[:, band] = formed_mass * np.exp(
            (1 - tx) * (1 - ty) * f00 + tx * (1 - ty) * f10
            + (1 - tx) * ty * f01 + tx * ty * f11)
    return out


def actual_particle_flux(data, observer_z):
    grid = renderer.get_fsps_grid(observer_z)
    flux = sb.star_flux(data["birth_tau"], data["metal"], data["initial_msun"],
                        data["a_snap"], data["a_grid"], data["tau_grid"], data["t_cum"],
                        grid[0], grid[1], grid[2])
    return flux * renderer.distance_factor(observer_z), grid


def selection_masks(data, observer_z, flux):
    """Recreate the paper renderer's six particle selections by projection."""
    n = int(round(2.0 * sb.STAMP_AS / sb.PIX_AS))
    xyz_pkpc = data["delta"] * data["a_snap"] / H * 1000.0
    kernels = [sb.psf_kernel(sb.PIX_AS, sb.FWHM_AS, sb.MOFFAT_BETA,
                             wing, sb.KERNEL_AS)[0] for wing in sb.WING_A]
    if abs(observer_z - data["z_phys"]) <= 1.0e-10:
        dm = float(renderer.COSMO_TNG.comoving_transverse_distance(observer_z).value)
        scale = 206265.0 / (H * dm)
    else:
        da = float(renderer.COSMO_TNG.angular_diameter_distance(observer_z).value)
        scale = 206265.0 * data["a_snap"] / (H * da)
    # Projection-specific axes are taken from the established renderer.
    output = []
    for axis_a, axis_b, los in LOS_AXES:
        ix = np.rint(data["delta"][:, axis_a] * scale / sb.PIX_AS).astype(np.int32) + n // 2
        iy = np.rint(data["delta"][:, axis_b] * scale / sb.PIX_AS).astype(np.int32) + n // 2
        on_image = (ix >= 0) & (ix < n) & (iy >= 0) & (iy < n)
        r2 = np.hypot(xyz_pkpc[:, axis_a], xyz_pkpc[:, axis_b])
        for iband, band in enumerate(BANDS):
            raw_img = sb.bin_flux(ix, iy, flux[:, iband], n)
            image = fftconvolve(raw_img, kernels[iband], mode="same")
            image[image < 0.0] = 0.0
            region, _ = sb.central_isophote(image, sb.THRESHOLDS[iband], sb.PIX_AS)
            iso = np.zeros(ix.size, dtype=bool)
            iso[on_image] = region[iy[on_image], ix[on_image]]
            masks = {
                "isophotal_particle_selection": iso,
                "b010_linked_group": data["link_masks"][0.10],
                "b012_linked_group": data["link_masks"][0.12],
                "b020_linked_group": data["link_masks"][0.20],
                "projected_aperture_30_pkpc": r2 <= 30.0,
                "Darwin_calibrated_sphere": data["radius3d_pkpc"] <= renderer.SPHERICAL_RADIUS_KPC[band],
            }
            output.append((los, band, iband, r2, masks))
    return output


def load_scores():
    scores = {}
    root = BASE / "growth_scores_r1audit_20261106"
    for rank in cohort_ranks():
        for row in read_csv(root / f"tree_rank_{rank:02d}_growth_scores.csv"):
            key = (int(row["sample_rank"]), int(row["snapshot"]), row["band"], row["method"],
                   row["line_of_sight"], round(float(row["observer_redshift"]), 6))
            if key in scores:
                raise RuntimeError(f"Duplicate existing score row {key}")
            scores[key] = row
    return scores


def class_from_start(start_raw, current_raw, start_snap):
    labels = np.full(current_raw["ids"].size, "not_central", dtype=object)
    current_c = np.flatnonzero(current_raw["central_valid"])
    idx, present = id_lookup(start_raw["ids"], current_raw["ids"][current_c])
    start_c_valid = np.zeros(current_c.size, dtype=bool)
    start_c_valid[present] = start_raw["central_valid"][idx[present]]
    persistent = start_c_valid
    new_by_formation = current_raw["birth_a"][current_c] > start_raw["a"] + 1.0e-7
    formed = ~persistent & new_by_formation
    old = ~persistent & ~new_by_formation
    labels[current_c[persistent]] = "persistent_since_snapshot50"
    labels[current_c[formed]] = "formed_since_snapshot50"
    old_ids = current_raw["ids"][current_c[old]]
    old_idx, old_present = id_lookup(start_raw["ids"], old_ids)
    old_labels = np.full(old_ids.size, "acquired_old_absent_from_parent_fof", dtype=object)
    old_labels[old_present & start_raw["satellite"][np.maximum(old_idx, 0)]] = (
        "acquired_old_from_resolved_satellite")
    non_sat = old_present & ~start_raw["satellite"][np.maximum(old_idx, 0)]
    old_labels[non_sat & ~start_raw["central"][np.maximum(old_idx, 0)]] = (
        "acquired_old_parent_fof_noncentral")
    old_labels[non_sat & start_raw["central"][np.maximum(old_idx, 0)]
               & ~start_raw["valid"][np.maximum(old_idx, 0)]] = (
        "acquired_old_invalid_start_membership")
    # If a central flag exists but was valid at the start, the ID would have
    # entered the persistent class above. Keep unexpected cases explicit.
    old_labels[non_sat & start_raw["central"][np.maximum(old_idx, 0)]
               & start_raw["valid"][np.maximum(old_idx, 0)]] = (
        "acquired_old_membership_audit_case")
    labels[current_c[old]] = old_labels
    return labels


def log_mass_ratio(a: float, b: float):
    if a <= 0 or b <= 0 or not np.isfinite(a + b):
        return float("nan")
    return float(np.log10(a / b))


def compute_size_metrics(rank, data, raw, start_raw, start_snap, actual_flux, actual_z,
                         score_index):
    central = data["central"]
    age = data["age_gyr"]
    z = np.clip(data["metal"], 1.0e-5, 0.05)
    mass = data["mass_msun"]
    formed = data["initial_msun"]
    central_w = mass[central]
    age_global = 10.0 ** mass_weighted_mean(np.log10(age[central]), central_w)
    z_global = 10.0 ** mass_weighted_mean(np.log10(z[central]), central_w)
    labels_raw = class_from_start(start_raw, raw, start_snap)
    class_use = labels_raw[data["use_idx"]]
    rows = []
    profile_rows = []
    for condition, observer_z in (("physical_observer_redshift", actual_z),
                                  ("fixed_observer_z0p197", OBSERVER_FIXED)):
        if condition == "physical_observer_redshift":
            full_flux, grid = actual_flux, renderer.get_fsps_grid(observer_z)
        else:
            full_flux, grid = actual_particle_flux(data, observer_z)
        no_age = particle_flux(np.full(age.size, age_global), z, formed, observer_z, grid)
        no_z = particle_flux(age, np.full(z.size, z_global), formed, observer_z, grid)
        no_both = particle_flux(np.full(age.size, age_global), np.full(z.size, z_global),
                                formed, observer_z, grid)
        flux_sets = {
            "full_stellar_ages_and_metallicities": full_flux,
            "uniform_age_varying_metallicity": no_age,
            "varying_age_uniform_metallicity": no_z,
            "uniform_age_and_metallicity": no_both,
        }
        for band_idx, band in enumerate(BANDS):
            for axis_a, axis_b, los in LOS_AXES:
                rproj = np.hypot(data["delta"][:, axis_a], data["delta"][:, axis_b])
                rproj *= data["a_snap"] / H * 1000.0
                r50_mass = weighted_r50(rproj, mass, central)
                r50_initial = weighted_r50(rproj, formed, central)
                r50_by_control = {
                    "current_stellar_mass": r50_mass,
                    "formed_stellar_mass": r50_initial,
                }
                for name, weights in flux_sets.items():
                    r50_by_control[name] = weighted_r50(rproj, weights[:, band_idx], central)
                row = {
                    "sample_rank": rank, "snapshot": data["snapshot"],
                    "physical_redshift": data["z_phys"], "observer_condition": condition,
                    "observer_redshift": observer_z, "band": band, "line_of_sight": los,
                    "central_mass_msun": float(mass[central].sum()),
                    "central_initial_mass_msun": float(formed[central].sum()),
                    "central_mass_weighted_mean_age_gyr": age_global,
                    "central_mass_weighted_geometric_mean_metallicity": z_global,
                    "n_central_stars": int(central.sum()),
                    **{f"R50_{name}_pkpc": value for name, value in r50_by_control.items()},
                }
                rows.append(row)

                # Compare full-SPS radii with the already released current
                # renderer for the matching condition and record discrepancies.
                if condition == "physical_observer_redshift":
                    expected = score_index.get((
                        rank, data["snapshot"], band, "central_bound_stars_reference", los,
                        round(observer_z, 6)))
                    if expected is None:
                        raise RuntimeError("Missing existing central-reference score")
                    expected_light = float(expected["central_bound_R50_pkpc"])
                    expected_mass = float(expected["central_bound_R50_mass_pkpc"])
                    row["existing_renderer_delta_logR50_light_dex"] = log_mass_ratio(
                        r50_by_control["full_stellar_ages_and_metallicities"], expected_light)
                    row["existing_renderer_delta_logR50_mass_dex"] = log_mass_ratio(
                        r50_mass, expected_mass)

                # Projected profiles are normalized by the current-mass R50.
                rel = rproj / r50_mass if r50_mass > 0 else np.full(rproj.size, np.nan)
                logz = np.log10(z)
                for ir, (lo, hi) in enumerate(zip(RADIAL_EDGES[:-1], RADIAL_EDGES[1:])):
                    annulus = central & (rel >= lo) & (rel < hi)
                    if not np.any(annulus):
                        continue
                    f = full_flux[:, band_idx]
                    m_sum, f_sum = float(mass[annulus].sum()), float(f[annulus].sum())
                    profile_rows.append({
                        "sample_rank": rank, "snapshot": data["snapshot"],
                        "physical_redshift": data["z_phys"], "observer_condition": condition,
                        "observer_redshift": observer_z, "band": band, "line_of_sight": los,
                        "radial_bin_R50mass": RADIAL_LABELS[ir],
                        "radial_lower_R50mass": lo,
                        "radial_upper_R50mass": "inf" if not np.isfinite(hi) else hi,
                        "n_stars": int(annulus.sum()), "current_stellar_mass_msun": m_sum,
                        "mass_weighted_mean_age_gyr": mass_weighted_mean(age[annulus], mass[annulus]),
                        "mass_weighted_mean_log10_metal_mass_fraction": mass_weighted_mean(
                            logz[annulus], mass[annulus]),
                        "log10_current_mass_to_light_ratio": (
                            math.log10(m_sum / f_sum) if m_sum > 0 and f_sum > 0 else ""),
                        "fraction_metallicity_clipped_low": float(np.mean(data["metal"][annulus] < 1e-5)),
                    })
    return rows, profile_rows


def summarize(values, rng):
    x = np.asarray(values, dtype=float)
    x = x[np.isfinite(x)]
    if not x.size:
        return {"n_branches": 0, "median": "", "p16": "", "p84": "",
                "median_ci95_low": "", "median_ci95_high": ""}
    draws = x[rng.integers(0, x.size, (BOOTSTRAPS, x.size))]
    ci = np.quantile(np.median(draws, axis=1), [0.025, 0.975])
    return {"n_branches": int(x.size), "median": float(np.median(x)),
            "p16": float(np.quantile(x, 0.16)), "p84": float(np.quantile(x, 0.84)),
            "median_ci95_low": float(ci[0]), "median_ci95_high": float(ci[1])}


def aggregate_budget(budget_rows, rng):
    quantity_map = {
        "net_central_mass_growth_over_start": "delta_over_start",
        "persistent_star_mass_change_over_start": "persistent_change_over_start",
        "preexisting_acquired_star_mass_over_start": "acquired_old_over_start",
        "newly_formed_star_mass_over_start": "new_formed_over_start",
        "central_stars_no_longer_bound_mass_over_start": "lost_over_start",
    }
    out = []
    groups = defaultdict(list)
    for row in budget_rows:
        groups[(int(row["snapshot_start"]), int(row["snapshot_end"]))].append(row)
    for (s0, s1), rows in sorted(groups.items()):
        for quantity, key in quantity_map.items():
            stats = summarize([float(r[key]) for r in rows], rng)
            out.append({"snapshot_start": s0, "snapshot_end": s1,
                        "quantity": quantity, "unit": "fraction of start central-subhalo stellar mass",
                        **stats})
        for quantity, key in (("closure_residual_fraction_of_start", "closure_residual_fraction_of_start"),):
            vals = np.abs([float(r[key]) for r in rows])
            out.append({"snapshot_start": s0, "snapshot_end": s1,
                        "quantity": quantity, "unit": "fraction of start mass",
                        "n_branches": len(vals), "median": float(np.median(vals)),
                        "p16": float(np.quantile(vals, .16)), "p84": float(np.quantile(vals, .84)),
                        "median_ci95_low": "", "median_ci95_high": ""})
    return out


def summarize_selection(method_rows, population_rows, rng):
    result = []
    for field in ("central_mass_recovery_fraction", "central_mass_fraction_of_selected",
                  "selected_total_mass_over_central_mass"):
        groups = defaultdict(lambda: defaultdict(list))
        for row in method_rows:
            key = (int(row["snapshot"]), row["band"], row["method"], "all_central_stars")
            groups[key][int(row["sample_rank"])].append(float(row[field]))
        for key, branches in groups.items():
            branch_values = [float(np.median(vals)) for vals in branches.values()]
            stats = summarize(branch_values, rng)
            result.append({"snapshot": key[0], "band": key[1], "method": key[2],
                           "population_class": key[3], "quantity": field, **stats})
    field = "capture_fraction_of_class"
    groups = defaultdict(lambda: defaultdict(list))
    for row in population_rows:
        key = (int(row["snapshot"]), row["band"], row["method"], row["population_class"])
        groups[key][int(row["sample_rank"])].append(float(row[field]))
    for key, branches in groups.items():
        branch_values = [float(np.median(vals)) for vals in branches.values()]
        stats = summarize(branch_values, rng)
        result.append({"snapshot": key[0], "band": key[1], "method": key[2],
                       "population_class": key[3], "quantity": field, **stats})
    return result


def build_growth_rows(size_rows):
    by_view = {(int(r["sample_rank"]), int(r["snapshot"]), r["band"],
                r["observer_condition"], r["line_of_sight"]): r for r in size_rows}
    out = []
    metrics = ("current_stellar_mass", "formed_stellar_mass",
               "full_stellar_ages_and_metallicities", "uniform_age_varying_metallicity",
               "varying_age_uniform_metallicity", "uniform_age_and_metallicity")
    for rank in cohort_ranks():
        for band in BANDS:
            for condition in ("physical_observer_redshift", "fixed_observer_z0p197"):
                per_metric_los = defaultdict(list)
                shapley_los = defaultdict(list)
                for los in ("x", "y", "z"):
                    a = by_view[(rank, 50, band, condition, los)]
                    b = by_view[(rank, 84, band, condition, los)]
                    for metric in metrics:
                        k = f"R50_{metric}_pkpc"
                        per_metric_los[metric].append(log_mass_ratio(float(b[k]), float(a[k])))
                    def lg(row, metric):
                        return math.log10(float(row[f"R50_{metric}_pkpc"]))
                    age_start = 0.5 * ((lg(a, "full_stellar_ages_and_metallicities")
                                        - lg(a, "uniform_age_varying_metallicity"))
                                       + (lg(a, "varying_age_uniform_metallicity")
                                          - lg(a, "uniform_age_and_metallicity")))
                    age_end = 0.5 * ((lg(b, "full_stellar_ages_and_metallicities")
                                      - lg(b, "uniform_age_varying_metallicity"))
                                     + (lg(b, "varying_age_uniform_metallicity")
                                        - lg(b, "uniform_age_and_metallicity")))
                    metal_start = 0.5 * ((lg(a, "full_stellar_ages_and_metallicities")
                                          - lg(a, "varying_age_uniform_metallicity"))
                                         + (lg(a, "uniform_age_varying_metallicity")
                                            - lg(a, "uniform_age_and_metallicity")))
                    metal_end = 0.5 * ((lg(b, "full_stellar_ages_and_metallicities")
                                        - lg(b, "varying_age_uniform_metallicity"))
                                       + (lg(b, "uniform_age_varying_metallicity")
                                          - lg(b, "uniform_age_and_metallicity")))
                    no_grad_growth = per_metric_los["uniform_age_and_metallicity"][-1]
                    full_growth = per_metric_los["full_stellar_ages_and_metallicities"][-1]
                    age_contribution = age_end - age_start
                    metal_contribution = metal_end - metal_start
                    shapley_los["age_gradient_contribution_dex"].append(age_contribution)
                    shapley_los["metallicity_gradient_contribution_dex"].append(metal_contribution)
                    shapley_los["uniform_population_growth_dex"].append(no_grad_growth)
                    shapley_los["decomposition_closure_dex"].append(
                        full_growth - (no_grad_growth + age_contribution + metal_contribution))
                for metric, vals in per_metric_los.items():
                    out.append({"sample_rank": rank, "band": band,
                                "observer_condition": condition,
                                "quantity": f"{metric}_growth_dex",
                                "branch_median_over_views": float(np.nanmedian(vals))})
                for metric, vals in shapley_los.items():
                    out.append({"sample_rank": rank, "band": band,
                                "observer_condition": condition, "quantity": metric,
                                "branch_median_over_views": float(np.nanmedian(vals))})
    return out


def summarize_growth(growth_rows, rng):
    groups = defaultdict(list)
    for row in growth_rows:
        groups[(row["band"], row["observer_condition"], row["quantity"])].append(
            float(row["branch_median_over_views"]))
    out = []
    for (band, condition, quantity), vals in sorted(groups.items()):
        out.append({"band": band, "observer_condition": condition,
                    "quantity": quantity, "unit": "dex", **summarize(vals, rng)})
    return out


def summarize_profiles(profile_rows, rng):
    # Reduce projections to one value per branch before any sample summary.
    grouped_view = defaultdict(list)
    fields = ("mass_weighted_mean_age_gyr", "mass_weighted_mean_log10_metal_mass_fraction",
              "log10_current_mass_to_light_ratio")
    for row in profile_rows:
        for field in fields:
            if row[field] == "":
                continue
            key = (int(row["sample_rank"]), int(row["snapshot"]), row["observer_condition"],
                   row["band"], row["radial_bin_R50mass"], field)
            grouped_view[key].append(float(row[field]))
    branch = defaultdict(list)
    for (rank, snap, condition, band, radial_bin, field), vals in grouped_view.items():
        branch[(snap, condition, band, radial_bin, field)].append(float(np.median(vals)))
    out = []
    for key, vals in sorted(branch.items()):
        out.append({"snapshot": key[0], "observer_condition": key[1], "band": key[2],
                    "radial_bin_R50mass": key[3], "quantity": key[4], **summarize(vals, rng)})
    return out


def make_figures(budget_summary, growth_summary):
    # Assembly budget figure: branch-level medians and branch bootstrap intervals.
    intervals = ((50, 67), (67, 84), (50, 84))
    quantities = (
        ("persistent_star_mass_change_over_start", "Persistent-star mass change"),
        ("preexisting_acquired_star_mass_over_start", "Acquired pre-existing stars"),
        ("newly_formed_star_mass_over_start", "Stars formed after start"),
        ("central_stars_no_longer_bound_mass_over_start", "No longer central"),
        ("net_central_mass_growth_over_start", "Net central mass growth"),
    )
    lookup = {(int(r["snapshot_start"]), int(r["snapshot_end"]), r["quantity"]): r
              for r in budget_summary}
    fig, ax = plt.subplots(figsize=(8.5, 4.8), constrained_layout=True)
    x = np.arange(len(intervals))
    width = 0.15
    for j, (quantity, label) in enumerate(quantities):
        med, lo, hi = [], [], []
        for s0, s1 in intervals:
            r = lookup[(s0, s1, quantity)]
            med.append(float(r["median"]))
            lo.append(float(r["median_ci95_low"]))
            hi.append(float(r["median_ci95_high"]))
        xpos = x + (j - (len(quantities) - 1) / 2) * width
        ax.bar(xpos, med, width=width, label=label, alpha=0.88)
        ax.errorbar(xpos, med, yerr=[np.array(med) - lo, hi - np.array(med)],
                    fmt="none", ecolor="black", elinewidth=0.7, capsize=2)
    ax.axhline(0, color="0.25", lw=0.8)
    ax.set_xticks(x, [r"$z\sim1.0$ to $0.5$", r"$z\sim0.5$ to $0.2$", r"$z\sim1.0$ to $0.2$"])
    ax.set_ylabel("Mass change / central stellar mass at start")
    ax.legend(ncol=2, fontsize=8, frameon=False)
    fig.savefig(OUTPUT_DIR / "assembly_mass_budget.png", dpi=200)
    plt.close(fig)

    # Controlled R50 growth figure. Each panel uses the same y range.
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.2), sharey=True, constrained_layout=True)
    categories = (
        ("current_stellar_mass_growth_dex", "Current-mass $R_{50}$"),
        ("full_stellar_ages_and_metallicities_growth_dex", "Full stellar populations"),
        ("uniform_age_varying_metallicity_growth_dex", "Uniform age; varying metallicity"),
        ("varying_age_uniform_metallicity_growth_dex", "Uniform metallicity; varying age"),
        ("uniform_age_and_metallicity_growth_dex", "Uniform age and metallicity"),
    )
    lookup = {(r["band"], r["observer_condition"], r["quantity"]): r
              for r in growth_summary}
    colors = plt.get_cmap("tab10").colors
    xpos = np.arange(2)
    width = 0.14
    for ib, band in enumerate(BANDS):
        ax = axes[ib]
        for j, (quantity, label) in enumerate(categories):
            med, lo, hi = [], [], []
            for condition in ("physical_observer_redshift", "fixed_observer_z0p197"):
                row = lookup[(band, condition, quantity)]
                med.append(float(row["median"]))
                lo.append(float(row["median_ci95_low"]))
                hi.append(float(row["median_ci95_high"]))
            xx = xpos + (j - (len(categories) - 1) / 2) * width
            ax.bar(xx, med, width=width, color=colors[j], label=label)
            ax.errorbar(xx, med, yerr=[np.array(med) - lo, hi - np.array(med)],
                        fmt="none", ecolor="black", elinewidth=0.65, capsize=2)
        ax.set_xticks(xpos, ["Physical observer\nredshifts", "Fixed observer\n$z=0.197$" ])
        ax.set_title(f"{band} band", loc="left", fontsize=10)
        ax.axhline(0, color="0.25", lw=0.8)
        ax.set_ylim(-0.2, 1.05)  # Reserve clear in-panel space for the legend.
    axes[0].set_ylabel(r"Median growth of projected $R_{50}$ [dex]")
    legend = axes[1].legend(fontsize=9, frameon=True, framealpha=0.94,
                            loc="upper right", borderaxespad=0.6)
    legend.set_in_layout(False)
    fig.savefig(OUTPUT_DIR / "size_growth_population_controls.png", dpi=200, bbox_inches="tight")
    plt.close(fig)


def main():
    global OUTPUT_DIR
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=OUT / "results")
    args = parser.parse_args()
    OUTPUT_DIR = args.output.resolve()
    if OUTPUT_DIR.exists() and any(OUTPUT_DIR.iterdir()):
        raise SystemExit(f"Refusing to overwrite existing analysis outputs: {OUTPUT_DIR}")
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    ranks = cohort_ranks()
    epochs = epoch_index()
    if len(ranks) != 36:
        raise RuntimeError(f"Expected 36 central-at-all-epochs branches, found {len(ranks)}")
    score_index = load_scores()
    kernels = None
    budget_rows, integrity_rows = [], []
    selected_rows, method_rows = [], []
    size_rows, profile_rows = [], []
    assembly_property_rows = []
    random = np.random.default_rng(SEED)

    for seq, rank in enumerate(ranks, 1):
        raw_by_snap, data_by_snap = {}, {}
        for snap in SNAPS:
            raw, data = load_epoch(rank, snap, epochs[rank][snap])
            raw_by_snap[snap], data_by_snap[snap] = raw, data
        # The exact central-particle identity ledger for adjacent intervals and
        # the complete z~1 to z~0.2 interval.
        transition_by_snap = {}
        for s0, s1 in ENDPOINT_INTERVALS:
            budget, transition = assembly_transition(
                raw_by_snap[s0], raw_by_snap[s1], s0, s1, rank)
            budget_rows.append(budget)
            if s0 == 50:
                transition_by_snap[s1] = transition
            if abs(float(budget["closure_residual_fraction_of_start"])) > 2e-12:
                raise RuntimeError(f"Mass ledger does not close for rank={rank}, {s0}->{s1}")

        # Audit ID continuity, initial-mass/birth-time stability, and central
        # particle masses against the SubLink catalogue at all outputs.
        for snap in SNAPS:
            raw, data = raw_by_snap[snap], data_by_snap[snap]
            integrity_rows.append({
                "sample_rank": rank, "snapshot": snap, "redshift": raw["z"],
                "n_parent_fof_stars": int(raw["ids"].size),
                "n_duplicate_particle_ids": int(raw["ids"].size - np.unique(raw["ids"]).size),
                "n_invalid_star_formation_times": int((~raw["valid"]).sum()),
                "n_central_members": int(raw["central_valid"].sum()),
                "n_central_type4_members_excluded_by_stellar_time_cut": int(
                    (raw["central"] & ~raw["valid"]).sum()),
                "central_particle_catalogue_mass_ratio": data["catalog_mass_ratio"],
                "n_central_satellite_mask_overlap": int((raw["central_valid"] & raw["satellite"]).sum()),
            })

        # Compute selection masks and validate the reproduced mass/count/R50
        # measurements against the frozen score table for every LOS and band.
        for snap in SNAPS:
            raw, data = raw_by_snap[snap], data_by_snap[snap]
            zobs = data["z_phys"]
            flux, grid = actual_particle_flux(data, zobs)
            labels_raw = (np.full(raw["ids"].size, "not_central", dtype=object)
                          if snap == 50 else class_from_start(
                              raw_by_snap[50], raw, 50))
            if snap == 50:
                labels_raw[raw["central_valid"]] = "baseline_central_at_snapshot50"
            class_use = labels_raw[data["use_idx"]]
            for los, band, iband, rproj, masks in selection_masks(data, zobs, flux):
                for method in METHODS:
                    chosen = masks[method]
                    central_chosen = chosen & data["central"]
                    total_mass = float(data["mass_msun"][chosen].sum())
                    central_mass = float(data["mass_msun"][data["central"]].sum())
                    selected_central_mass = float(data["mass_msun"][central_chosen].sum())
                    r50 = weighted_r50(rproj, flux[:, iband], chosen)
                    expected = score_index[(rank, snap, band, method, los, round(zobs, 6))]
                    score_count = int(float(expected["selected_particle_count"]))
                    score_mass = float(expected["stellar_mass_msun"])
                    score_r50 = float(expected["projected_R50_pkpc"])
                    delta_logmass = log_mass_ratio(total_mass, score_mass)
                    delta_logr = log_mass_ratio(r50, score_r50)
                    if (int(chosen.sum()) != score_count or abs(delta_logmass) > 2e-7
                            or abs(delta_logr) > 2e-7):
                        raise RuntimeError(
                            f"Selection reproduction mismatch rank={rank} snap={snap} "
                            f"{band}/{method}/{los}: count {chosen.sum()} vs {score_count}; "
                            f"dlogM={delta_logmass}; dlogR={delta_logr}")
                    method_rows.append({
                        "sample_rank": rank, "snapshot": snap, "physical_redshift": data["z_phys"],
                        "band": band, "line_of_sight": los, "method": method,
                        "selected_particle_count": int(chosen.sum()),
                        "selected_total_mass_msun": total_mass,
                        "selected_central_mass_msun": selected_central_mass,
                        "central_mass_msun": central_mass,
                        "central_mass_recovery_fraction": selected_central_mass / central_mass,
                        "central_mass_fraction_of_selected": (
                            selected_central_mass / total_mass if total_mass > 0 else ""),
                        "selected_total_mass_over_central_mass": total_mass / central_mass,
                        "projected_half_light_radius_pkpc": r50,
                        "delta_log10_mass_vs_existing_score": delta_logmass,
                        "delta_log10_R50_vs_existing_score": delta_logr,
                    })
                    for label in np.unique(class_use[data["central"]]):
                        pop = data["central"] & (class_use == label)
                        pop_mass = float(data["mass_msun"][pop].sum())
                        if pop_mass <= 0:
                            continue
                        selected_pop_mass = float(data["mass_msun"][pop & chosen].sum())
                        selected_rows.append({
                            "sample_rank": rank, "snapshot": snap,
                            "physical_redshift": data["z_phys"], "band": band,
                            "line_of_sight": los, "method": method,
                            "population_class": str(label), "population_mass_msun": pop_mass,
                            "selected_population_mass_msun": selected_pop_mass,
                            "capture_fraction_of_class": selected_pop_mass / pop_mass,
                            "selected_population_fraction_of_central": selected_pop_mass / central_mass,
                            "population_particle_count": int(pop.sum()),
                            "selected_population_particle_count": int((pop & chosen).sum()),
                        })
            # Size analysis uses exact central membership and reports radial
            # age/Z/M/L profiles plus age/metallicity-gradient controls.
            start = raw_by_snap[50]
            size_part, profile_part = compute_size_metrics(
                rank, data, raw, start, 50, flux, zobs, score_index)
            size_rows.extend(size_part)
            profile_rows.extend(profile_part)
        print(f"rank {rank} of 36 complete", flush=True)

    # Reproduction error summary for the existing six selections.
    max_mass_delta = max(abs(float(r["delta_log10_mass_vs_existing_score"]))
                         for r in method_rows)
    max_r50_delta = max(abs(float(r["delta_log10_R50_vs_existing_score"]))
                        for r in method_rows)
    if max_mass_delta > 2e-7 or max_r50_delta > 2e-7:
        raise RuntimeError("Reproduced selections differ from the frozen measurements")

    budget_summary = aggregate_budget(budget_rows, random)
    selection_summary = summarize_selection(method_rows, selected_rows, random)
    growth_rows = build_growth_rows(size_rows)
    growth_summary = summarize_growth(growth_rows, random)
    profile_summary = summarize_profiles(profile_rows, random)

    write_csv(OUTPUT_DIR / "assembly_mass_budget_per_branch.csv", budget_rows)
    write_csv(OUTPUT_DIR / "assembly_mass_budget_summary.csv", budget_summary)
    write_csv(OUTPUT_DIR / "input_membership_integrity.csv", integrity_rows)
    write_csv(OUTPUT_DIR / "selection_metrics_per_view.csv", method_rows)
    write_csv(OUTPUT_DIR / "selection_population_capture_per_view.csv", selected_rows)
    write_csv(OUTPUT_DIR / "selection_population_capture_summary.csv", selection_summary)
    write_csv(OUTPUT_DIR / "central_size_metrics_per_view.csv", size_rows)
    write_csv(OUTPUT_DIR / "central_size_growth_per_branch.csv", growth_rows)
    write_csv(OUTPUT_DIR / "central_size_growth_summary.csv", growth_summary)
    write_csv(OUTPUT_DIR / "radial_population_profiles_per_view.csv", profile_rows)
    write_csv(OUTPUT_DIR / "radial_population_profiles_summary.csv", profile_summary)
    make_figures(budget_summary, growth_summary)

    ledger_abs = [abs(float(r["closure_residual_fraction_of_start"])) for r in budget_rows]
    score_delta = {"max_abs_log_mass_difference_dex": max_mass_delta,
                   "max_abs_log_R50_difference_dex": max_r50_delta,
                   "selection_rows_compared": len(method_rows)}
    manifest = {
        "simulation": "TNG300-1", "sample": "36 central-at-all-targets SubLink branches",
        "snapshots": {"50": 0.9972942257819404, "67": 0.503047523244883,
                      "84": 0.19728418237600986},
        "input_directory": "../tree_epoch_membership",
        "input_hdf5_count": 108,
        "input_fields": ["Coordinates", "Masses", "ParticleIDs", "GFM_InitialMass",
                         "GFM_Metallicity", "GFM_StellarFormationTime",
                         "CentralSubhaloParticle", "RemoveSatelliteGalaxy"],
        "mass_unit_conversion": "10^10 Msun/h with h=0.6774",
        "assembly_ledger": (
            "central mass change = current-mass change of persistent central stars + "
            "current mass of pre-existing stars newly assigned to the central subhalo + "
            "current mass of stars formed since the start output and central at the end - "
            "start-output mass of stars no longer assigned to the central subhalo"),
        "formation_classification": (
            "GFM_StellarFormationTime scale factor after the start output marks stars formed "
            "between outputs; the particle birth host is not known, so the class is not labelled in situ"),
        "size_controls": {
            "full": "particle ages and metallicities retained",
            "uniform_age_varying_metallicity": "mass-weighted mean log age replaces particle ages",
            "varying_age_uniform_metallicity": "mass-weighted mean log metallicity replaces particle metallicities",
            "uniform_age_and_metallicity": "both population variables held spatially uniform",
            "observer_conditions": ["physical observer redshift", "fixed observer z=0.197284"],
            "radius_definition": "projected R50 using unblurred particle fluxes for exact central-subhalo members",
        },
        "statistical_unit": "descendant branch; three projections reduced to branch median",
        "bootstrap_replicates": BOOTSTRAPS, "bootstrap_seed": SEED,
        "maximum_mass_ledger_closure_fraction": max(ledger_abs),
        "existing_renderer_reproduction": score_delta,
        "data_limit": (
            "The cutouts contain no per-particle velocity or gravitational-potential field. "
            "Satellite mask is boolean and does not identify a donor subhalo."),
        "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    }
    write_json(OUTPUT_DIR / "analysis_manifest.json", manifest)
    readme = f"""# TNG300 assembly and population-weighting analysis

This analysis follows the locked 36-branch TNG300-1 sample at snapshots 50, 67,
and 84. The analysis uses the parent-FoF stellar cutouts already included in
`../tree_epoch_membership/`. No new simulation calculation or network fetch was
required.

## Analysis 1: particle-ID mass ledger

For each interval, the ledger matches `ParticleIDs` and decomposes the central
Subfind stellar-mass change into the current-mass change of stars assigned to
the central subhalo at both endpoints, pre-existing stars assigned to the
central subhalo only at the later endpoint, stars formed between endpoints
and assigned to the later central subhalo, and the earlier mass of stars no
longer assigned to the central subhalo. The ledger closes to a maximum absolute
fractional residual of `{max(ledger_abs):.3e}` of the starting central mass.

The earlier `RemoveSatelliteGalaxy` flag identifies stars in resolved
non-central subhaloes with at least 100 stellar particles. The flag is boolean,
so a source subhalo ID is unavailable. Pre-existing stars missing from the
earlier parent-FoF cutout have unknown prior locations. The analysis does not
label stars formed between outputs as in situ because the birth host is not
provided. The current cutouts contain no particle velocities or potential.

The six selection masks from the paper renderer were recreated at the physical
observer redshift and compared against the saved per-particle score tables.
Maximum absolute differences are `{max_mass_delta:.3e}` dex in selected mass and
`{max_r50_delta:.3e}` dex in projected half-light radius across
`{len(method_rows)}` branch, snapshot, band, projection, and method records.

## Analysis 2: stellar-population size controls

The analysis measures projected half-light radii for central-subhalo stars
using the existing dust-free SSP flux grids. The full calculation retains each
particle's age and total metal mass fraction. Three counterfactual controls
replace the particle ages, metallicities, or both with the central population's
mass-weighted mean log value at each output. The controls preserve particle
positions and formed masses. The Shapley split of the log-radius difference
shares the age-metallicity interaction equally between the two terms.

The radii use full particle fluxes at particle positions, as in the released
central-reference measurement. The calculation excludes image PSF and
isophotal selection and therefore isolates population weighting. The physical
observer-redshift condition changes the observer-frame filter between outputs.
The fixed-redshift condition keeps the observer redshift at 0.197284.

Radial profiles use projected radius in units of the current-mass half-mass
radius. Each projection is reduced to a branch median before sample summaries.
Intervals in the summary tables are 95 per cent percentile bootstrap intervals
over the 36 descendant branches, using {BOOTSTRAPS} resamples and seed {SEED}.

## Outputs

- `assembly_mass_budget_per_branch.csv` and `assembly_mass_budget_summary.csv`
- `selection_metrics_per_view.csv`
- `selection_population_capture_per_view.csv` and
  `selection_population_capture_summary.csv`
- `central_size_metrics_per_view.csv`, `central_size_growth_per_branch.csv`,
  and `central_size_growth_summary.csv`
- `radial_population_profiles_per_view.csv` and
  `radial_population_profiles_summary.csv`
- `assembly_mass_budget.png` and `size_growth_population_controls.png`
- `analysis_manifest.json`

The TNG data documentation defines `GFM_StellarFormationTime` as the stellar
formation scale factor, `GFM_InitialMass` as birth mass, and
`GFM_Metallicity` as total metal mass fraction. See the
[TNG data specifications](https://www.illustris-project.org/data/docs/specifications/).
"""
    with (OUTPUT_DIR / "README.md").open("x") as stream:
        stream.write(readme)

    print(json.dumps({
        "output_directory": str(OUTPUT_DIR),
        "n_branches": len(ranks), "n_budget_rows": len(budget_rows),
        "n_selection_rows": len(method_rows), "n_size_view_rows": len(size_rows),
        "max_ledger_closure_fraction": max(ledger_abs),
        "selection_reproduction": score_delta,
        "completed": True,
    }, indent=2), flush=True)


if __name__ == "__main__":
    main()
