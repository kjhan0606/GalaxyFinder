#!/usr/bin/env python3
"""Counterfactual redshift experiment for one fixed Darwin-1 halo.

The z=1.99 particle distribution and stellar ages are frozen. The same halo is
rendered at archived observed redshifts using four controlled photometric
normalizations. This isolates rendering effects; it is not a model of the
halo's physical evolution or passive stellar aging.
"""
from __future__ import annotations

import csv
import argparse
import json
import sys
from pathlib import Path

import numpy as np
from scipy.signal import fftconvolve

_BASE = Path(__file__).resolve().parents[2]
ROOT = _BASE.parent if _BASE.name == "PGalF_major_revision_20261013" else _BASE
PROJECT = ROOT / "PGalF_major_revision_20261013"
WORK = PROJECT / "artificial_redshift_20261029"
SOURCE = PROJECT / "scripts" / "source"
SUPP = PROJECT / "inputs" / "supplement"
sys.path.insert(0, str(SOURCE))
import sb_vs_b as sb  # noqa: E402

OBS_GRIDS = {
    "0.190": SUPP / "fsps_z0.190.npz",
    "0.227": SUPP / "fsps_z0.227.npz",
    "0.500": SUPP / "fsps_z0.500.npz",
    "1.002": SUPP / "fsps_z1.002.npz",
    "1.990": SUPP / "fsps_z1.990.npz",
}
REST_GRID = SUPP / "fsps_rest_z0.npz"
B_TEST = (0.05, 0.10, 0.20)


def read_stars_photometry_only(f, counts):
    """Read image fields without decoding particle IDs.

    The legacy 00967 stellar records have ambiguous ID32/ID64 padding. This
    experiment does not use particle IDs; positions, masses, ages, and
    metallicities occupy the same audited offsets in both record layouts.
    The archived-IoU gate below checks that the resulting image calculation
    reproduces the original measurements.
    """
    f.seek(int(counts[0]) * sb.DM + int(counts[1]) * sb.GAS
           + int(counts[2]) * sb.SINK, 1)
    n = int(counts[3])
    if n == 0:
        return None
    raw = f.read(n * sb.STAR)
    if len(raw) != n * sb.STAR:
        raise ValueError("Short stellar block in photometry-only reader")
    offsets = {"x": 0, "y": 8, "z": 16, "mass": 48,
               "tp": 72, "zp": 80, "mass0": 112}
    out = {
        key: np.ndarray(n, dtype="<f8", buffer=raw, offset=offset,
                        strides=(sb.STAR,)).copy()
        for key, offset in offsets.items()
    }
    if not all(np.all(np.isfinite(values)) for values in out.values()):
        raise ValueError("Non-finite image field in stellar records")
    return out


def load_grid(path: Path):
    with np.load(path) as data:
        return float(data["z"]), data["age"].copy(), data["logz"].copy(), data["flux"].copy()


def unwrap(q, centre, box):
    delta = q - centre
    return delta - box * np.round(delta / box)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--step", required=True)
    parser.add_argument("--nstar", required=True, type=int)
    parser.add_argument("--tag", required=True,
                        help="Unique result suffix, e.g. static_population_redshift_v4")
    args = parser.parse_args()
    step = str(args.step).zfill(5)
    nstar = args.nstar
    result = WORK / "results" / f"{step}_{nstar}_{args.tag}"
    if result.exists():
        if any(result.iterdir()):
            raise FileExistsError(f"Refusing to overwrite nonempty {result}")
        # A failed read can leave only the empty directory created below; it
        # contains no analysis output and can be safely reused on a retry.
        result.rmdir()
    result.mkdir(parents=True)

    if step not in sb.STEPS:
        raise RuntimeError(f"Unknown snapshot step {step}")
    a_sim, gpath, bpath, wanted = sb.STEPS[step]
    cache_tag = step
    if step == "00967":
        # The locked first-paper sample at this output was measured from the
        # stellar-link-off top100/0 catalogue, not top100/1 used by the later
        # satellite study. Use its matching background list and offset cache.
        gpath = PROJECT.parent / "Darwin1/top100/0/GALFIND.DATA.00967"
        bpath = PROJECT.parent / "Darwin1/top100/0/background_ptl.00967"
        cache_tag = "00967_off"
    halos = sb.cached_offsets(gpath, False, cache_tag)
    locked = {
        tuple(line.split())
        for line in (PROJECT / "inputs" / "halos.txt").read_text().splitlines()
        if line.strip()
    }
    if (step, str(nstar)) not in locked:
        raise RuntimeError(f"Halo {step}:{nstar} is not in the locked revision sample")
    entry = next((row for row in halos if row[1] is not None and int(row[1]) == nstar), None)
    if entry is None:
        raise RuntimeError(f"Halo {step}:{nstar} is not in the archived sample")
    ih, nstar, _mstar, halo_offset = entry
    if step == "00967":
        metadata_path = PROJECT / "inputs" / "supplement" / "halo_metadata.tsv"
        metadata_lines = [line for line in metadata_path.read_text().splitlines()
                          if line and not line.startswith("#")]
        metadata_rows = list(csv.DictReader(metadata_lines, delimiter="\t"))
        metadata_row = next(r for r in metadata_rows
                            if r["step"] == step and int(r["nstar"]) == nstar)
        if (int(metadata_row["catalogue_row"]) != int(ih)
                or int(metadata_row["byte_offset"]) != int(halo_offset)):
            raise RuntimeError("00967 input does not match the locked /0 halo metadata")
    backs = sb.cached_offsets(bpath, True, cache_tag)
    if ih >= len(backs):
        raise RuntimeError(f"No matching background block for catalogue row {ih}")

    if step == "00967":
        # IDs are neither read nor used in this image-only experiment. The
        # old 00967 writer leaves enough padding ambiguity to defeat a safe
        # ID-width inference, so decode only the independently fixed fields.
        sb.read_stars = read_stars_photometry_only
    pieces = sb.load_halo(gpath, int(halo_offset))
    background = sb.load_halo(bpath, int(backs[ih][3]), background=True)[0]
    cat, centre, n_sat, n_satellite_galaxies = sb.parent_stars(pieces, background)
    box = 44.171267
    coords = tuple(unwrap(cat[k], centre[j], box) for j, k in enumerate(("x", "y", "z")))
    mass = cat["mass"]
    if step == "00967" and len(mass) != int(metadata_row["n_parent"]):
        raise RuntimeError(
            f"00967 parent-star count differs from locked metadata: "
            f"{len(mass)} != {metadata_row['n_parent']}"
        )
    formed_mass = np.clip(cat["mass0"], 0.0, None) / (sb.H0 / 100.0)
    a_grid, tau_grid, t_cum = sb.build_time(a_sim)
    rest_z, age_grid, logz_grid, rest_flux_grid = load_grid(REST_GRID)
    if rest_z != 0.0:
        raise RuntimeError("Unexpected redshift in rest-frame SSP grid")

    membership_path = PROJECT / "results" / f"{step}_{nstar}" / "membership.npz"
    keep_by_b = {}
    with np.load(membership_path) as archived_membership:
        if int(archived_membership["n_stars"]) != len(mass):
            raise RuntimeError("Archived FoF membership length does not match reconstructed parent stars")
        for b in B_TEST:
            key = f"b{b:.2f}"
            keep = np.unpackbits(archived_membership[key], count=len(mass)).astype(bool)
            keep_by_b[b] = keep
            print("REUSED_ARCHIVED_FOF", b, int(keep.sum()), flush=True)

    audit = json.loads((PROJECT / "results" / "photometry_audit.json").read_text())
    rows: list[dict[str, object]] = []
    baseline_ref = {
        (r["los"], r["band"], round(float(r["b"]), 2)): float(r["iou"])
        for r in csv.DictReader((PROJECT / "results" / "recomputed_curves.tsv").open(), delimiter="\t")
        if r["step"] == step and int(r["nstar"]) == nstar
        and r["variant"] == "consistent_distance"
    }

    kernels = [
        sb.psf_kernel(sb.PIX_AS, sb.FWHM_AS, sb.MOFFAT_BETA, amp, sb.KERNEL_AS)[0]
        for amp in sb.WING_A
    ]
    n_pix = int(round(2.0 * sb.STAMP_AS / sb.PIX_AS))
    h = sb.H0 / 100.0
    for grid_name, grid_path in OBS_GRIDS.items():
        z_obs, z_age, z_logz, observer_flux_grid = load_grid(grid_path)
        # Keep the stellar population fixed at its z=1.99 ages. Only the
        # apparent spectrum, cosmological dimming, and angular scale change.
        observer_flux = sb.star_flux(
            cat["tp"], cat["zp"], formed_mass, a_sim, a_grid, tau_grid, t_cum,
            z_age, z_logz, observer_flux_grid,
        )
        rest_flux = sb.star_flux(
            cat["tp"], cat["zp"], formed_mass, a_sim, a_grid, tau_grid, t_cum,
            age_grid, logz_grid, rest_flux_grid,
        )
        audit_row = next(r for r in audit["records"] if abs(float(r["z"]) - z_obs) < 1e-8)
        dl_pc = float(audit_row["fsps_distance_pc"])
        d_sim_pc = float(audit_row["simulation_distance_pc"])
        consistent_distance_factor = (dl_pc / d_sim_pc) ** 2
        distance_only_scale = (10.0 / d_sim_pc) ** 2
        observer_flux *= consistent_distance_factor
        modes = {
            "observer_full": observer_flux,
            "observer_no_tolman_factor": observer_flux * (1.0 + z_obs) ** 4,
            "rest_sed_with_distance": rest_flux * distance_only_scale,
            "rest_sed_no_dimming": rest_flux,
        }

        # Input positions are comoving Mpc/h at z_sim. Preserve their
        # physical size at z_sim when projecting at z_obs: r_proper=a_sim*r_com,
        # and D_A=D_M/(1+z_obs). This recovers the archived scale at z_obs=z_sim.
        dm = sb.COSMO.comoving_transverse_distance(z_obs).to_value("Mpc")
        scale = a_sim * (1.0 + z_obs) * 206265.0 / (h * dm)
        pix_coords = []
        for j, k in sb.LOS_AXES:
            ix = np.rint(coords[j] * scale / sb.PIX_AS).astype(np.int32) + n_pix // 2
            iy = np.rint(coords[k] * scale / sb.PIX_AS).astype(np.int32) + n_pix // 2
            on = (ix >= 0) & (ix < n_pix) & (iy >= 0) & (iy < n_pix)
            pix_coords.append((ix, iy, on))

        for mode_name, weights in modes.items():
            for ib, band in enumerate(sb.BAND_NAME):
                for ia, (los, (ix, iy, on)) in enumerate(zip(sb.LOS_NAME, pix_coords, strict=True)):
                    image = sb.bin_flux(ix, iy, weights[:, ib], n_pix)
                    image = fftconvolve(image, kernels[ib], mode="same")
                    image[image < 0.0] = 0.0
                    region, area = sb.central_isophote(image, sb.THRESHOLDS[ib], sb.PIX_AS)
                    inside = np.zeros(len(weights), dtype=bool)
                    inside[on] = region[iy[on], ix[on]]
                    for b, keep in keep_by_b.items():
                        w = weights[:, ib]
                        l_group = float(w[keep].sum())
                        l_inside = float(w[inside].sum())
                        l_both = float(w[keep & inside].sum())
                        denominator = l_group + l_inside - l_both
                        iou = l_both / denominator if denominator > 0 else float("nan")
                        completeness = l_both / l_inside if l_inside > 0 else float("nan")
                        precision = l_both / l_group if l_group > 0 else float("nan")
                        rows.append({
            "step": step, "nstar": nstar, "z_sim": 1.0 / a_sim - 1.0,
                            "z_obs": z_obs, "grid": grid_name, "flux_mode": mode_name,
                            "los": los, "band": band, "b": b, "n_parent": len(mass),
                            "n_satellite_stars_masked": n_sat,
                            "n_satellite_galaxies_masked": n_satellite_galaxies,
                            "n_target_stars": int(inside.sum()), "area_arcsec2": area,
                            "R_eq_arcsec": float(np.sqrt(area / np.pi)) if area > 0 else 0.0,
                            "L_group": l_group, "L_inside": l_inside, "L_common": l_both,
                            "iou": iou, "completeness": completeness, "precision": precision,
                        })
                    del image, region, inside
            del weights
        del observer_flux, rest_flux
        print("RENDERED", grid_name, "z_obs", z_obs, flush=True)

    # The z_obs=z_sim observer-frame rendering must reproduce the archived
    # boundary scores before the controlled comparisons are accepted.
    checks = []
    z_sim = 1.0 / a_sim - 1.0
    for row in rows:
        if row["flux_mode"] != "observer_full" or abs(float(row["z_obs"]) - z_sim) > 1e-5:
            continue
        reference = baseline_ref[(str(row["los"]), str(row["band"]), round(float(row["b"]), 2))]
        delta = float(row["iou"]) - reference
        checks.append(delta)
        if abs(delta) > 5.1e-5:
            raise RuntimeError(f"Baseline reproduction failed: {row}; archived IoU={reference}")
    if len(checks) != 18:
        raise RuntimeError(f"Expected 18 baseline checks, got {len(checks)}")

    out = result / "metrics.tsv"
    with out.open("x", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    metadata = {
        "experiment": "fixed z=1.99 Darwin-1 halo, static stellar ages and particle distribution",
        "step": step, "nstar": nstar, "z_sim": z_sim,
        "source_catalogue": str(gpath), "source_background": str(bpath),
        "locked_sample_row": ({"catalogue_row": int(metadata_row["catalogue_row"]),
                                "byte_offset": int(metadata_row["byte_offset"]),
                                "n_parent": int(metadata_row["n_parent"])}
                               if step == "00967" else None),
        "observation_grids": {k: str(v) for k, v in OBS_GRIDS.items()},
        "photometric_modes": {
            "observer_full": "redshifted FSPS spectrum with FSPS-to-Darwin luminosity-distance correction",
            "observer_no_tolman_factor": "observer_full multiplied by (1+z)^4; retains redshifted bandpass",
            "rest_sed_with_distance": "z=0 FSPS spectrum scaled by the Darwin luminosity distance; removes bandpass shift",
            "rest_sed_no_dimming": "z=0 absolute FSPS spectrum without distance scaling",
        },
        "distance_normalization": "all observer-frame and distance-scaled rest-spectrum images use the Darwin luminosity distance; observer spectra retain their FSPS/Darwin correction factor",
        "geometry": {"coordinate_convention": "proper size fixed at z_sim; angle = a_sim * r_comoving / D_A(z_obs)",
                     "angular_scale_factor_relative_to_archived": "a_sim * (1 + z_obs)"},
        "fixed_instrument": {"pixel_arcsec": sb.PIX_AS, "PSF_FWHM_arcsec": sb.FWHM_AS,
                             "thresholds_mag_arcsec2": list(sb.THRESHOLDS)},
        "baseline_reproduction_n": len(checks),
        "baseline_max_abs_iou_difference": max(abs(x) for x in checks),
        "stellar_field_reader": ("photometry-only; IDs intentionally not decoded"
                                 if step == "00967" else "verified standard reader"),
        "limits": [
            "This is a counterfactual rendering experiment, not halo evolution.",
            "Stellar ages, particle masses, and proper positions are frozen at z_sim; comoving coordinates are converted with a_sim before projection.",
            "The fixed angular instrument means proper pixel and PSF scales vary with z_obs.",
            "No sky noise, dust, IGM absorption, or passive stellar aging is included.",
        ],
    }
    (result / "metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")
    print(json.dumps(metadata, indent=2), flush=True)


if __name__ == "__main__":
    main()
