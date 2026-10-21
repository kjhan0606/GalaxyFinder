#!/usr/bin/env python3
"""Stream stellar metallicity and kinematic summaries from NewGalFinder output.

The 112-byte HaloInfo and SubInfo records below match NewGalFinder/defs.h;
the particle records match NewGalFinder/ramses.h and the build flags in
NewGalFinder/Makefile (dptype=double, idtype=int, NCHEM=9, NDUST=4, no
OUTPUT_PARTICLE_POTENTIAL, no NBODY), including C alignment padding. The
NewDD reader scales particle velocities to km/s before writing these records;
the extractor preserves those values without a second conversion. The stellar-
particle field zp stores the total birth-metallicity mass fraction. The program
reports number-weighted and current-mass-weighted means in that convention; it
does not convert them to solar units or to [Fe/H]. The extractor also measures
each galaxy's current stellar half-mass radius about its stellar mass centroid.
For duplicate audits, it uses abs(ID) as the tracking key documented by NewDD.
"""

from __future__ import annotations

import argparse
import csv
import math
import os
import re
from pathlib import Path

import numpy as np


HALO_DTYPE = np.dtype(
    {
        "names": [
            "nsub", "ndm", "nstar", "nsink", "ngas", "npall",
            "totm", "mdm", "mgas", "msink", "mstar",
            "x", "y", "z", "vx", "vy", "vz",
        ],
        "formats": ["<i4"] * 6 + ["<f8"] * 11,
        "offsets": list(range(0, 24, 4)) + list(range(24, 112, 8)),
        "itemsize": 112,
    }
)

SUB_DTYPE = np.dtype(
    {
        "names": [
            "npdm", "npgas", "npsink", "npstar", "npall",
            "totm", "mdm", "mgas", "msink", "mstar",
            "x", "y", "z", "vx", "vy", "vz",
        ],
        "formats": ["<i4"] * 5 + ["<f8"] * 11,
        "offsets": list(range(0, 20, 4)) + list(range(24, 112, 8)),
        "itemsize": 112,
    }
)

# PmType/StarType, from ramses.h.  Six phase-space coordinates, current mass,
# integer ID and particle tags, birth time, total metallicity zp, nine element
# abundances, and initial mass.  Native C padding is represented explicitly.
STAR_DTYPE = np.dtype(
    {
        "names": ["x", "y", "z", "vx", "vy", "vz", "mass", "id",
                  "levelp", "family", "tag", "tp", "zp", "chem", "mass0"],
        # The six bytes between tag and tp are native C alignment padding,
        # already represented by the explicit field offsets below.  They
        # must not appear as an extra dtype field here.
        "formats": ["<f8"] * 7 + ["<i4", "<i4", "i1", "i1",
                                   "<f8", "<f8", ("<f8", (9,)), "<f8"],
        "offsets": [0, 8, 16, 24, 32, 40, 48, 56, 60, 64, 65, 72, 80, 88, 160],
        "itemsize": 168,
    }
)

# These two records are skipped, but their exact sizes are needed to locate
# the following stellar records in each subhalo payload.
DM_DTYPE = STAR_DTYPE
SINK_DTYPE = np.dtype({"names": ["bytes"], "formats": ["V168"], "itemsize": 168})
GAS_DTYPE = np.dtype({"names": ["bytes"], "formats": ["V152"], "itemsize": 152})


def read_record(stream, dtype: np.dtype, label: str):
    raw = stream.read(dtype.itemsize)
    if not raw:
        return None
    if len(raw) != dtype.itemsize:
        raise EOFError(f"truncated {label} at byte {stream.tell() - len(raw)}")
    return np.frombuffer(raw, dtype=dtype, count=1)[0]


def read_star_chunk(stream, count: int) -> np.ndarray:
    raw = stream.read(count * STAR_DTYPE.itemsize)
    expected = count * STAR_DTYPE.itemsize
    if len(raw) != expected:
        raise EOFError(f"truncated star payload: read {len(raw)} of {expected} bytes")
    return np.frombuffer(raw, dtype=STAR_DTYPE, count=count)


def particle_tracking_ids(ids: np.ndarray) -> np.ndarray:
    """Return the identity key used downstream for signed NewDD particle IDs."""
    # NewDD/rd_part.c preserves the signed RAMSES ID but documents abs(ID) as
    # the downstream tracking key. Widen before abs() to avoid int32 overflow
    # for the valid ID -2**31.
    return np.abs(np.asarray(ids, dtype=np.int64))


def read_cosmology_values(
    input_path: Path, snapshot: str
) -> tuple[float, float, float, float, float, float, float]:
    info_path = (input_path.parents[2] / f"output_{snapshot}" /
                 f"info_{snapshot}.txt")
    if not info_path.is_file():
        raise FileNotFoundError(f"missing RAMSES info file {info_path}")
    values = {}
    with info_path.open(encoding="ascii") as stream:
        for line in stream:
            match = re.match(
                r"\s*(aexp|H0|omega_m|omega_l|omega_k|unit_l)\s*=\s*([0-9.eE+-]+)",
                line,
            )
            if match:
                values[match.group(1)] = float(match.group(2))
    required = {"aexp", "H0", "omega_m", "omega_l", "unit_l"}
    if not required.issubset(values):
        raise ValueError(f"cosmology fields missing from {info_path}: {sorted(required - values.keys())}")
    omega_k = values.get("omega_k", 1.0 - values["omega_m"] - values["omega_l"])
    hubble_h = values["H0"] / 100.0
    if values["aexp"] <= 0 or hubble_h <= 0 or values["unit_l"] <= 0:
        raise ValueError(f"non-positive cosmology scale in {info_path}")
    # RAMSES unit_l is the proper length of one code unit. In a cosmological
    # run, unit_l = a * Lbox(cMpc/h) * Mpc / h.
    boxlen_cmpc_h = values["unit_l"] * hubble_h / (values["aexp"] * 3.085678e24)
    if not math.isfinite(boxlen_cmpc_h) or boxlen_cmpc_h <= 0:
        raise ValueError(f"invalid comoving box length in {info_path}")
    return (
        values["aexp"], hubble_h, boxlen_cmpc_h, values["H0"],
        values["omega_m"], values["omega_l"], omega_k,
    )


class RamsesTimeConverter:
    """Convert RAMSES cosmological stellar birth coordinates to proper age.

    RAMSES `part2birth` maps the stored stellar `tp` coordinate to proper time
    by interpolating its Friedmann tables.  This builds a denser lookup table
    from the same expansion-history equations so the conversion can be applied
    to particle chunks without treating `tp` as either a scale factor or an
    elapsed age.
    """

    def __init__(self, aexp, H0, omega_m, omega_l, omega_k, ntable=100_001):
        amin = 1.0e-3
        if not (0 < aexp <= 1.0) or H0 <= 0 or omega_m <= 0:
            raise ValueError("invalid cosmology for stellar birth-time conversion")
        if ntable < 1001:
            raise ValueError("stellar-age lookup table is too coarse")

        scale = np.geomspace(amin, 1.0, ntable, dtype=np.float64)
        log_scale = np.log(scale)
        radicand = omega_m + omega_l * scale**3 + omega_k * scale
        dadt = np.sqrt(radicand / scale)
        dadtau = np.sqrt(scale**3 * radicand)
        if not (np.isfinite(dadt).all() and np.isfinite(dadtau).all()
                and (dadt > 0).all() and (dadtau > 0).all()):
            raise ValueError("RAMSES Friedmann expansion rates are invalid")

        proper_rate = scale / dadt
        birth_rate = scale / dadtau
        dlog_scale = np.diff(log_scale)
        proper_from_amin = np.zeros_like(scale)
        birth_from_amin = np.zeros_like(scale)
        proper_from_amin[1:] = np.cumsum(
            0.5 * (proper_rate[:-1] + proper_rate[1:]) * dlog_scale,
            dtype=np.float64,
        )
        birth_from_amin[1:] = np.cumsum(
            0.5 * (birth_rate[:-1] + birth_rate[1:]) * dlog_scale,
            dtype=np.float64,
        )

        # RAMSES stores both time coordinates as look-back values relative to
        # a=1; they therefore increase from negative values to zero in `scale`.
        self.tau = birth_from_amin - birth_from_amin[-1]
        self.proper_time = proper_from_amin - proper_from_amin[-1]
        self.scale = scale
        self.aexp = float(aexp)
        self.tau_at_output = float(np.interp(aexp, scale, self.tau))
        self.proper_time_at_output = float(np.interp(aexp, scale, self.proper_time))
        # Match the H0-to-seconds conversion used by RAMSES `part2birth`.
        self.hubble_time_gyr = (
            3.08e24 / (H0 * 1.0e5) / (365.0 * 24.0 * 3600.0 * 1.0e9)
        )

    def ages_from_birth_coordinate(self, tp: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        birth = np.asarray(tp, dtype=np.float64)
        valid = (
            np.isfinite(birth)
            & (birth >= self.tau[0])
            & (birth <= self.tau_at_output + 1.0e-8)
        )
        ages = np.full(birth.shape, np.nan, dtype=np.float64)
        if valid.any():
            birth_time = np.interp(
                np.minimum(birth[valid], self.tau_at_output), self.tau,
                self.proper_time,
            )
            values = (self.proper_time_at_output - birth_time) * self.hubble_time_gyr
            good_age = np.isfinite(values) & (values >= -1.0e-8)
            valid_indices = np.flatnonzero(valid)
            valid[valid_indices[~good_age]] = False
            ages[valid_indices[good_age]] = np.maximum(values[good_age], 0.0)
        return ages, valid


class StarAccumulator:
    """Stellar moments and a periodic, mass-weighted 3D size measurement."""

    def __init__(self, boxlen_cmpc_h, physical_kpc_per_cmpc_h, time_converter):
        self.pivot = None
        self.position_anchor = None
        self.boxlen_cmpc_h = float(boxlen_cmpc_h)
        self.physical_kpc_per_cmpc_h = float(physical_kpc_per_cmpc_h)
        self.time_converter = time_converter
        self.stellar_particle_current_mass_sum = 0.0
        self.age_weight = 0.0
        self.age_weighted_sum = 0.0
        self.nstar_age_valid = 0
        self.weight = 0.0
        self.weighted_offset = np.zeros(3, dtype=np.float64)
        self.weighted_offset2 = np.zeros(3, dtype=np.float64)
        self.number = 0
        self.number_mass = 0
        self.number_offset = np.zeros(3, dtype=np.float64)
        self.number_offset2 = np.zeros(3, dtype=np.float64)
        self.z_weight = 0.0
        self.z_weighted_sum = 0.0
        self.z_count = 0
        self.z_number_sum = 0.0
        self.z_number_count = 0
        self.star_ids = []
        self.stellar_offsets = []
        self.stellar_position_masses = []
        self.nstar_position_mass_valid = 0
        self.kinematic_offsets = []
        self.kinematic_velocities = []
        self.kinematic_masses = []

    def add(self, stars: np.ndarray):
        if not len(stars):
            return
        self.star_ids.append(stars["id"].copy())
        velocity = np.column_stack((stars["vx"], stars["vy"], stars["vz"]))
        mass = stars["mass"].astype(np.float64, copy=False)
        position = np.column_stack((stars["x"], stars["y"], stars["z"]))
        valid_current_mass = np.isfinite(mass) & (mass > 0)
        if valid_current_mass.any():
            self.stellar_particle_current_mass_sum += float(
                mass[valid_current_mass].sum(dtype=np.float64)
            )
        age, valid_birth_time = self.time_converter.ages_from_birth_coordinate(stars["tp"])
        valid_age = valid_current_mass & valid_birth_time
        if valid_age.any():
            age_mass = mass[valid_age]
            self.age_weight += float(age_mass.sum(dtype=np.float64))
            self.age_weighted_sum += float(
                (age_mass * age[valid_age]).sum(dtype=np.float64)
            )
            self.nstar_age_valid += int(valid_age.sum())
        valid_position_mass = (np.isfinite(position).all(axis=1)
                               & np.isfinite(mass) & (mass > 0))
        if valid_position_mass.any():
            self.nstar_position_mass_valid += int(valid_position_mass.sum())
            halfbox = 0.5 * self.boxlen_cmpc_h
            valid_positions = position[valid_position_mass]
            if self.position_anchor is None:
                self.position_anchor = np.mod(valid_positions[0], self.boxlen_cmpc_h)
            offset = (valid_positions - self.position_anchor + halfbox) % self.boxlen_cmpc_h - halfbox
            self.stellar_offsets.append(offset)
            self.stellar_position_masses.append(mass[valid_position_mass].copy())
            valid_velocity = np.isfinite(velocity).all(axis=1)[valid_position_mass]
            if valid_velocity.any():
                self.kinematic_offsets.append(offset[valid_velocity])
                self.kinematic_velocities.append(velocity[valid_position_mass][valid_velocity].copy())
                self.kinematic_masses.append(mass[valid_position_mass][valid_velocity].copy())
        z = stars["zp"].astype(np.float64, copy=False)
        valid_z_number = np.isfinite(z)
        if valid_z_number.any():
            self.z_number_sum += float(z[valid_z_number].sum(dtype=np.float64))
            self.z_number_count += int(valid_z_number.sum())
        valid_z_mass = np.isfinite(mass) & (mass > 0) & np.isfinite(z)
        if valid_z_mass.any():
            wz = mass[valid_z_mass]
            self.z_weight += float(wz.sum(dtype=np.float64))
            self.z_weighted_sum += float((wz * z[valid_z_mass]).sum(dtype=np.float64))
            self.z_count += int(valid_z_mass.sum())
        valid_v = np.isfinite(velocity).all(axis=1)
        valid_w = valid_v & np.isfinite(mass) & (mass > 0)
        if not valid_v.any():
            return
        if self.pivot is None:
            self.pivot = velocity[np.flatnonzero(valid_v)[0]].copy()

        v = velocity[valid_v]
        off = v - self.pivot
        self.number += len(v)
        self.number_offset += off.sum(axis=0, dtype=np.float64)
        self.number_offset2 += np.square(off).sum(axis=0, dtype=np.float64)

        if valid_w.any():
            vw = velocity[valid_w]
            w = mass[valid_w]
            self.number_mass += int(valid_w.sum())
            offw = vw - self.pivot
            self.weight += float(w.sum(dtype=np.float64))
            self.weighted_offset += (offw * w[:, None]).sum(axis=0, dtype=np.float64)
            self.weighted_offset2 += (np.square(offw) * w[:, None]).sum(axis=0, dtype=np.float64)

    @staticmethod
    def _sigmas(sum1, sum2, weight):
        if weight <= 0:
            return (math.nan,) * 5
        mean = sum1 / weight
        variance = np.maximum(sum2 / weight - np.square(mean), 0.0)
        components = np.sqrt(variance)
        return (float(np.sqrt(variance.mean())), float(np.linalg.norm(components)),
                float(components[0]), float(components[1]), float(components[2]))

    def result(self):
        n_id = sum(len(x) for x in self.star_ids)
        n_id_unique = (
            int(np.unique(particle_tracking_ids(np.concatenate(self.star_ids))).size)
            if self.star_ids else 0
        )
        age_result = {
            "nstar_age_valid": self.nstar_age_valid,
            "mean_stellar_age_current_mass_weighted_gyr": (
                self.age_weighted_sum / self.age_weight if self.age_weight > 0 else math.nan
            ),
        }
        if self.number == 0 or self.pivot is None:
            return {
                "nstar_velocity_valid": self.number,
                "nstar_velocity_mass_valid": self.number_mass,
                "nstar_position_mass_valid": self.nstar_position_mass_valid,
                "nstar_velocity_within_rhalf": 0,
                "nstar_id_count": n_id,
                "nstar_id_unique": n_id_unique,
                "nstar_id_duplicate": n_id - n_id_unique,
                "stellar_particle_current_mass_sum_msun_h": self.stellar_particle_current_mass_sum,
                **age_result,
                "nstar_birth_metallicity_valid": self.z_count,
                "mean_birth_metallicity_number_weighted": math.nan,
                "mstar_velocity_weight": self.weight,
                "mean_birth_metallicity_current_mass_weighted": math.nan,
                "mean_vx_star_kms": math.nan,
                "mean_vy_star_kms": math.nan,
                "mean_vz_star_kms": math.nan,
                "sigma_star_1d_mass_weighted_kms": math.nan,
                "sigma_star_3d_mass_weighted_kms": math.nan,
                "sigma_star_x_mass_weighted_kms": math.nan,
                "sigma_star_y_mass_weighted_kms": math.nan,
                "sigma_star_z_mass_weighted_kms": math.nan,
                "sigma_star_1d_number_weighted_kms": math.nan,
                "sigma_star_3d_number_weighted_kms": math.nan,
                "sigma_star_x_number_weighted_kms": math.nan,
                "sigma_star_y_number_weighted_kms": math.nan,
                "sigma_star_z_number_weighted_kms": math.nan,
                "rhalf3d_stellar_mass_weighted_kpc": math.nan,
                "sigma_star_1d_mass_weighted_within_rhalf_kms": math.nan,
            }

        weighted = self._sigmas(self.weighted_offset, self.weighted_offset2, self.weight)
        numbered = self._sigmas(self.number_offset, self.number_offset2, self.number)
        vmean = self.pivot + self.weighted_offset / self.weight if self.weight else [math.nan] * 3
        rhalf_kpc = sigma_rhalf = math.nan
        nstar_velocity_within_rhalf = 0
        if self.stellar_offsets:
            offsets = np.concatenate(self.stellar_offsets)
            position_mass = np.concatenate(self.stellar_position_masses)
            total_position_mass = float(position_mass.sum(dtype=np.float64))
            if total_position_mass > 0:
                stellar_centre = np.average(offsets, axis=0, weights=position_mass)
                halfbox = 0.5 * self.boxlen_cmpc_h
                centered = (offsets - stellar_centre + halfbox) % self.boxlen_cmpc_h - halfbox
                radius_cmpc_h = np.linalg.norm(centered, axis=1)
                order = np.argsort(radius_cmpc_h, kind="mergesort")
                cumulative_mass = np.cumsum(position_mass[order], dtype=np.float64)
                radius_index = int(np.searchsorted(cumulative_mass,
                                                   0.5 * total_position_mass,
                                                   side="left"))
                rhalf_cmpc_h = float(radius_cmpc_h[order[min(radius_index, len(order) - 1)]])
                rhalf_kpc = rhalf_cmpc_h * self.physical_kpc_per_cmpc_h
                if self.kinematic_offsets:
                    kin_offsets = np.concatenate(self.kinematic_offsets)
                    kin_velocities = np.concatenate(self.kinematic_velocities)
                    kin_masses = np.concatenate(self.kinematic_masses)
                    centered_kin = (kin_offsets - stellar_centre + halfbox) % self.boxlen_cmpc_h - halfbox
                    inside = np.linalg.norm(centered_kin, axis=1) <= rhalf_cmpc_h
                    if inside.any():
                        nstar_velocity_within_rhalf = int(inside.sum())
                        aperture_mass = kin_masses[inside]
                        aperture_velocity = kin_velocities[inside]
                        aperture_mean = np.average(aperture_velocity, axis=0, weights=aperture_mass)
                        aperture_variance = np.average(
                            np.square(aperture_velocity - aperture_mean),
                            axis=0, weights=aperture_mass,
                        )
                        sigma_rhalf = float(np.sqrt(aperture_variance.mean()))
        return {
            "nstar_velocity_valid": self.number,
            "nstar_velocity_mass_valid": self.number_mass,
            "nstar_position_mass_valid": self.nstar_position_mass_valid,
            "nstar_velocity_within_rhalf": nstar_velocity_within_rhalf,
            "nstar_id_count": n_id,
            "nstar_id_unique": n_id_unique,
            "nstar_id_duplicate": n_id - n_id_unique,
            "stellar_particle_current_mass_sum_msun_h": self.stellar_particle_current_mass_sum,
            **age_result,
            "nstar_birth_metallicity_valid": self.z_count,
            "mean_birth_metallicity_number_weighted": self.z_number_sum / self.z_number_count
            if self.z_number_count else math.nan,
            "mstar_velocity_weight": self.weight,
            "mean_birth_metallicity_current_mass_weighted": self.z_weighted_sum / self.z_weight if self.z_weight else math.nan,
            "mean_vx_star_kms": float(vmean[0]),
            "mean_vy_star_kms": float(vmean[1]),
            "mean_vz_star_kms": float(vmean[2]),
            "sigma_star_1d_mass_weighted_kms": weighted[0],
            "sigma_star_3d_mass_weighted_kms": weighted[1],
            "sigma_star_x_mass_weighted_kms": weighted[2],
            "sigma_star_y_mass_weighted_kms": weighted[3],
            "sigma_star_z_mass_weighted_kms": weighted[4],
            "sigma_star_1d_number_weighted_kms": numbered[0],
            "sigma_star_3d_number_weighted_kms": numbered[1],
            "sigma_star_x_number_weighted_kms": numbered[2],
            "sigma_star_y_number_weighted_kms": numbered[3],
            "sigma_star_z_number_weighted_kms": numbered[4],
            "rhalf3d_stellar_mass_weighted_kpc": rhalf_kpc,
            "sigma_star_1d_mass_weighted_within_rhalf_kms": sigma_rhalf,
        }


def extract(input_path: Path, output_path: Path, chunk_size: int):
    if output_path.exists():
        raise FileExistsError(f"refusing to overwrite {output_path}")
    temp_output_path = output_path.with_name(
        f".{output_path.name}.partial.{os.getpid()}"
    )
    if temp_output_path.exists():
        raise FileExistsError(f"refusing to overwrite {temp_output_path}")
    input_size = input_path.stat().st_size
    match = re.search(r"GALFIND\.DATA\.(\d{5})$", input_path.name)
    snapshot = match.group(1) if match else "unknown"
    if snapshot == "unknown":
        raise ValueError(f"cannot infer snapshot number from {input_path.name}")
    (scale_factor, hubble_h, boxlen_cmpc_h, H0, omega_m, omega_l,
     omega_k) = read_cosmology_values(input_path, snapshot)
    time_converter = RamsesTimeConverter(
        scale_factor, H0, omega_m, omega_l, omega_k,
    )
    redshift = 1.0 / scale_factor - 1.0
    columns = [
        "snapshot", "scale_factor", "redshift", "hubble_h", "host_index", "sub_index",
        "is_central_max_stellar_mass", "host_nsub",
        "host_fof_dm_star_mass_msun_h", "host_stellar_mass_msun_h",
        "host_x_cmpc_h", "host_y_cmpc_h", "host_z_cmpc_h",
        "host_vx_kms", "host_vy_kms", "host_vz_kms",
        "galaxy_nstar", "galaxy_stellar_mass_msun_h",
        "stellar_particle_current_mass_sum_msun_h",
        "stellar_particle_mass_sum_over_catalogue_mstar",
        "galaxy_dm_star_mass_msun_h",
        "galaxy_dark_matter_mass_msun_h", "galaxy_gas_mass_msun_h",
        "galaxy_sink_mass_msun_h",
        "galaxy_x_cmpc_h", "galaxy_y_cmpc_h", "galaxy_z_cmpc_h",
        "galaxy_vx_kms", "galaxy_vy_kms", "galaxy_vz_kms",
        "nstar_velocity_valid", "nstar_velocity_mass_valid",
        "nstar_id_count", "nstar_id_unique",
        "nstar_position_mass_valid", "nstar_velocity_within_rhalf",
        "nstar_id_duplicate", "nstar_ids_in_multi_galaxy_groups",
        "host_duplicate_star_id_excess", "host_multi_galaxy_duplicate_id_count",
        "nstar_age_valid", "mean_stellar_age_current_mass_weighted_gyr",
        "nstar_birth_metallicity_valid", "mean_birth_metallicity_number_weighted",
        "mstar_velocity_weight", "mean_birth_metallicity_current_mass_weighted",
        "mean_vx_star_kms", "mean_vy_star_kms",
        "mean_vz_star_kms", "sigma_star_1d_mass_weighted_kms",
        "sigma_star_3d_mass_weighted_kms", "sigma_star_x_mass_weighted_kms",
        "sigma_star_y_mass_weighted_kms", "sigma_star_z_mass_weighted_kms",
        "sigma_star_1d_number_weighted_kms",
        "sigma_star_3d_number_weighted_kms", "sigma_star_x_number_weighted_kms",
        "sigma_star_y_number_weighted_kms", "sigma_star_z_number_weighted_kms",
        "rhalf3d_stellar_mass_weighted_kpc",
        "sigma_star_1d_mass_weighted_within_rhalf_kms",
    ]
    with input_path.open("rb") as src, temp_output_path.open("x", newline="") as dst:
        writer = csv.DictWriter(dst, fieldnames=columns)
        writer.writeheader()
        host_index = 0
        galaxy_count = 0
        # Keep the global cross-host ID audit bounded in memory. NewGalFinder
        # stores signed 32-bit IDs; NewDD documents abs(ID) as the downstream
        # tracking key. A bitset covers every possible absolute key, including
        # abs(-2**31), in about 256 MiB, independent of catalogue size.
        global_id_bitmap = np.zeros((1 << 28) + 1, dtype=np.uint8)
        global_unique_count = 0
        while True:
            host = read_record(src, HALO_DTYPE, "HaloInfo")
            if host is None:
                break
            nsub = int(host["nsub"])
            if nsub < 0 or nsub > 1_000_000:
                raise ValueError(f"implausible nsub={nsub} at byte {src.tell() - HALO_DTYPE.itemsize}")
            host_rows = []
            for sub_index in range(nsub):
                sub = read_record(src, SUB_DTYPE, "SubInfo")
                if sub is None:
                    raise EOFError(f"missing SubInfo for host {host_index}, sub {sub_index}")
                counts = {k: int(sub[k]) for k in ("npdm", "npgas", "npsink", "npstar", "npall")}
                if any(v < 0 for v in counts.values()):
                    raise ValueError(f"negative member count at host {host_index}, sub {sub_index}: {counts}")
                if sum(counts[k] for k in ("npdm", "npgas", "npsink", "npstar")) != counts["npall"]:
                    raise ValueError(f"component counts do not sum to npall at host {host_index}, sub {sub_index}: {counts}")
                src.seek(counts["npdm"] * DM_DTYPE.itemsize
                         + counts["npgas"] * GAS_DTYPE.itemsize
                         + counts["npsink"] * SINK_DTYPE.itemsize, 1)
                if src.tell() > input_size:
                    raise EOFError(f"nonstellar payload exceeds file at host {host_index}, sub {sub_index}")
                acc = StarAccumulator(
                    boxlen_cmpc_h, scale_factor / hubble_h * 1000.0,
                    time_converter,
                )
                remaining = counts["npstar"]
                while remaining:
                    nread = min(remaining, chunk_size)
                    acc.add(read_star_chunk(src, nread))
                    remaining -= nread
                row = {
                    "snapshot": snapshot,
                    "scale_factor": scale_factor,
                    "redshift": redshift,
                    "hubble_h": hubble_h,
                    "host_index": host_index,
                    "sub_index": sub_index,
                    "host_nsub": nsub,
                    "host_fof_dm_star_mass_msun_h": float(host["totm"]),
                    "host_stellar_mass_msun_h": float(host["mstar"]),
                    "host_x_cmpc_h": float(host["x"]),
                    "host_y_cmpc_h": float(host["y"]),
                    "host_z_cmpc_h": float(host["z"]),
                    "host_vx_kms": float(host["vx"]),
                    "host_vy_kms": float(host["vy"]),
                    "host_vz_kms": float(host["vz"]),
                    "galaxy_nstar": counts["npstar"],
                    "galaxy_stellar_mass_msun_h": float(sub["mstar"]),
                    "galaxy_dm_star_mass_msun_h": float(sub["totm"]),
                    "galaxy_dark_matter_mass_msun_h": float(sub["mdm"]),
                    "galaxy_gas_mass_msun_h": float(sub["mgas"]),
                    "galaxy_sink_mass_msun_h": float(sub["msink"]),
                    "galaxy_x_cmpc_h": float(sub["x"]),
                    "galaxy_y_cmpc_h": float(sub["y"]),
                    "galaxy_z_cmpc_h": float(sub["z"]),
                    "galaxy_vx_kms": float(sub["vx"]),
                    "galaxy_vy_kms": float(sub["vy"]),
                    "galaxy_vz_kms": float(sub["vz"]),
                }
                row.update(acc.result())
                catalogue_mstar = float(row["galaxy_stellar_mass_msun_h"])
                particle_mstar = float(row["stellar_particle_current_mass_sum_msun_h"])
                row["stellar_particle_mass_sum_over_catalogue_mstar"] = (
                    particle_mstar / catalogue_mstar if catalogue_mstar > 0 else math.nan
                )
                row["_star_ids"] = (np.concatenate(acc.star_ids)
                                    if acc.star_ids else np.empty(0, dtype=np.int32))
                host_rows.append(row)
            if host_rows:
                # Check whether a star ID is assigned to more than one galaxy
                # in this FoF host.  The per-galaxy accumulator already counts
                # repeated IDs within a single galaxy; this host-level pass
                # finds the separate cross-galaxy failure mode.
                id_arrays = [row["_star_ids"] for row in host_rows if row["_star_ids"].size]
                if id_arrays:
                    owner_arrays = [np.full(row["_star_ids"].size, i, dtype=np.int32)
                                    for i, row in enumerate(host_rows)
                                    if row["_star_ids"].size]
                    ids = particle_tracking_ids(np.concatenate(id_arrays))
                    # Tracking keys repeated within this host are handled by the
                    # per-galaxy and cross-galaxy checks below. Only unique
                    # host IDs enter the global bitset, which detects reuse in
                    # a different FoF host without retaining every ID until
                    # the end of the snapshot.
                    unique_ids, inverse, id_counts = np.unique(
                        ids, return_inverse=True, return_counts=True
                    )
                    unsigned_ids = unique_ids.astype(np.uint64, copy=False)
                    byte_index = (unsigned_ids >> 3).astype(np.intp)
                    bit_mask = np.left_shift(
                        np.uint8(1), (unsigned_ids & np.uint64(7)).astype(np.uint8)
                    )
                    seen_before = (global_id_bitmap[byte_index] & bit_mask) != 0
                    if seen_before.any():
                        sample_ids = unique_ids[seen_before][:20].tolist()
                        temp_output_path.unlink(missing_ok=True)
                        raise ValueError(
                            "stellar particle tracking keys (abs(ID)) occur in multiple FoF hosts: "
                            f"n_duplicate_ids={int(seen_before.sum())}, "
                            f"sample_ids={sample_ids}; refusing to publish catalogue"
                        )
                    np.bitwise_or.at(global_id_bitmap, byte_index, bit_mask)
                    global_unique_count += int(unique_ids.size)
                    owners = np.concatenate(owner_arrays)
                    host_excess = int(ids.size - unique_ids.size)
                    duplicate_ids = np.flatnonzero(id_counts > 1)
                    if duplicate_ids.size:
                        pairs, pair_counts = np.unique(
                            np.column_stack((inverse, owners)), axis=0, return_counts=True
                        )
                        pair_id = pairs[:, 0]
                        pair_owner = pairs[:, 1]
                        owner_counts = np.bincount(pair_id, minlength=len(unique_ids))
                        multi_owner_pairs = owner_counts[pair_id] > 1
                        multi_member_counts = np.bincount(
                            pair_owner[multi_owner_pairs],
                            weights=pair_counts[multi_owner_pairs],
                            minlength=len(host_rows),
                        ).astype(np.int64)
                        multi_id_count = int(np.count_nonzero(owner_counts > 1))
                    else:
                        multi_member_counts = np.zeros(len(host_rows), dtype=np.int64)
                        multi_id_count = 0
                else:
                    host_excess = 0
                    multi_id_count = 0
                    multi_member_counts = np.zeros(len(host_rows), dtype=np.int64)
                for i, row in enumerate(host_rows):
                    row["nstar_ids_in_multi_galaxy_groups"] = int(multi_member_counts[i])
                    row["host_duplicate_star_id_excess"] = host_excess
                    row["host_multi_galaxy_duplicate_id_count"] = multi_id_count
                    del row["_star_ids"]
                # Match NewGalFinder/src/binding.c: central_stellar_core,
                # which identifies the core with the largest stellar mass.
                central_idx = max(
                    range(len(host_rows)),
                    key=lambda i: host_rows[i]["galaxy_stellar_mass_msun_h"],
                )
                for i, row in enumerate(host_rows):
                    row["is_central_max_stellar_mass"] = int(i == central_idx)
                    writer.writerow(row)
                galaxy_count += len(host_rows)
            host_index += 1
    # Publish only a fully written catalogue. Linking fails atomically if a
    # final output appeared concurrently, so existing results are preserved.
    os.link(temp_output_path, output_path)
    temp_output_path.unlink()
    print(f"snapshot={snapshot} hosts={host_index} galaxies={galaxy_count} "
          f"global_unique_star_ids={global_unique_count} output={output_path}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path, help="GALFIND.DATA.NNNNN file")
    parser.add_argument("output", type=Path, help="new CSV output path; existing files are never overwritten")
    parser.add_argument("--chunk-size", type=int, default=250_000,
                        help="star records processed per chunk (default: 250000)")
    args = parser.parse_args()
    if args.chunk_size <= 0:
        parser.error("--chunk-size must be positive")
    if not args.input.is_file():
        parser.error(f"input does not exist: {args.input}")
    extract(args.input, args.output, args.chunk_size)


if __name__ == "__main__":
    main()
