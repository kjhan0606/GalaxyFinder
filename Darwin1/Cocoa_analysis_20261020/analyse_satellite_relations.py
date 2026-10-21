#!/usr/bin/env python3
"""Exploratory satellite metallicity/kinematics analysis for ASCII extracts.

The internal-dispersion model includes host fixed effects and stellar mass.
As a separate descriptive check, the host orbital-dispersion comparison pairs
the metal-poor and metal-rich satellites within each host by stellar mass. It
concerns galaxy-centre orbital velocities, not the primary internal stellar
velocity dispersion; missing centre velocities do not remove a satellite from
the internal-dispersion sample.
"""

from __future__ import annotations

import argparse
import csv
import math
import os
from pathlib import Path

import numpy as np


def rank_average(x: np.ndarray) -> np.ndarray:
    order = np.argsort(x, kind="mergesort")
    sorted_x = x[order]
    ranks = np.empty(len(x), dtype=np.float64)
    i = 0
    while i < len(x):
        j = i + 1
        while j < len(x) and sorted_x[j] == sorted_x[i]:
            j += 1
        ranks[order[i:j]] = 0.5 * (i + 1 + j)
        i = j
    return ranks


def corr(x, y):
    if len(x) < 3:
        return math.nan
    rx, ry = rank_average(x), rank_average(y)
    return float(np.corrcoef(rx, ry)[0, 1])


def within_host(values: np.ndarray, hosts: np.ndarray) -> np.ndarray:
    out = values.astype(np.float64, copy=True)
    _, inverse, counts = np.unique(hosts, return_inverse=True, return_counts=True)
    sums = np.bincount(inverse, weights=out)
    out -= sums[inverse] / counts[inverse]
    return out


def fixed_host_ols(xz, xm, y, hosts):
    X = np.column_stack((within_host(xz, hosts), within_host(xm, hosts)))
    yw = within_host(y, hosts)
    good = np.isfinite(X).all(axis=1) & np.isfinite(yw)
    X, yw, h = X[good], yw[good], hosts[good]
    n, k, ng = len(yw), X.shape[1], len(np.unique(h))
    if ng < 2 or n <= ng + k or np.linalg.matrix_rank(X) < k:
        return np.full(k, math.nan), np.full(k, math.nan), n, ng
    beta = np.linalg.lstsq(X, yw, rcond=None)[0]
    resid = yw - X @ beta
    inv = np.linalg.pinv(X.T @ X)
    meat = np.zeros((2, 2), dtype=np.float64)
    groups = np.unique(h)
    for g in groups:
        m = h == g
        score = X[m].T @ resid[m]
        meat += np.outer(score, score)
    # The finite-sample CR1 correction must count both slopes and the
    # host-specific intercepts in the residual degrees of freedom.
    correction = (ng / (ng - 1)) * ((n - 1) / (n - ng - k))
    cov = correction * inv @ meat @ inv
    return beta, np.sqrt(np.maximum(np.diag(cov), 0.0)), n, ng


def fixed_host_ols_controls(predictors, y, hosts):
    """Host fixed-effects OLS with cluster-robust errors and extra controls."""
    X = np.column_stack([within_host(x, hosts) for x in predictors])
    yw = within_host(y, hosts)
    good = np.isfinite(X).all(axis=1) & np.isfinite(yw)
    X, yw, h = X[good], yw[good], hosts[good]
    n, k, ng = len(yw), X.shape[1], len(np.unique(h))
    if (ng < 2 or n <= ng + k or np.linalg.matrix_rank(X) < k):
        return np.full(X.shape[1], math.nan), np.full(X.shape[1], math.nan), len(yw), len(np.unique(h))
    beta = np.linalg.lstsq(X, yw, rcond=None)[0]
    resid = yw - X @ beta
    bread = np.linalg.pinv(X.T @ X)
    meat = np.zeros((X.shape[1], X.shape[1]), dtype=np.float64)
    groups = np.unique(h)
    for g in groups:
        m = h == g
        score = X[m].T @ resid[m]
        meat += np.outer(score, score)
    correction = (ng / (ng - 1)) * ((n - 1) / (n - ng - k))
    cov = correction * bread @ meat @ bread
    return beta, np.sqrt(np.maximum(np.diag(cov), 0.0)), n, ng


def fixed_host_massbin_slope(xz, xm, y, hosts, width=0.3, min_cell_n=3):
    massbin = np.floor((xm - np.nanmin(xm)) / width).astype(np.int64)
    cells = np.rec.fromarrays((hosts, massbin), names=("host", "massbin"))
    _, inverse, counts = np.unique(cells, return_inverse=True, return_counts=True)
    keep = counts[inverse] >= min_cell_n
    x, yy, hh, cc = xz[keep], y[keep], hosts[keep], inverse[keep]
    if len(x) < 10:
        return math.nan, math.nan, int(len(x)), int(len(np.unique(hh)))
    _, ci, cn = np.unique(cc, return_inverse=True, return_counts=True)
    sx = np.bincount(ci, weights=x)
    sy = np.bincount(ci, weights=yy)
    xw = x - sx[ci] / cn[ci]
    yw = yy - sy[ci] / cn[ci]
    xx = float(xw @ xw)
    if xx <= 0:
        return math.nan, math.nan, int(len(x)), int(len(np.unique(hh)))
    n, ncell = len(xw), len(np.unique(cc))
    if n <= ncell + 1:
        return math.nan, math.nan, int(n), int(len(np.unique(hh)))
    beta = float(xw @ yw / xx)
    resid = yw - beta * xw
    host_ids, hi = np.unique(hh, return_inverse=True)
    scores = np.bincount(hi, weights=xw * resid)
    ng = len(host_ids)
    if ng < 2:
        return beta, math.nan, int(n), int(ng)
    correction = (ng / (ng - 1)) * ((n - 1) / (n - ncell - 1))
    se = math.sqrt(correction * float(scores @ scores) / xx**2)
    return beta, se, n, ng


def host_bootstrap_interval(xz, xm, y, hosts, nboot, rng, massbin_width=None,
                            min_cell_n=3):
    """Percentile interval from resampling whole host haloes."""
    host_ids, host_inverse = np.unique(hosts, return_inverse=True)
    nhost = len(host_ids)
    if nhost < 2 or nboot <= 0:
        return math.nan, math.nan

    if massbin_width is None:
        xz_w = within_host(xz, hosts)
        xm_w = within_host(xm, hosts)
        y_w = within_host(y, hosts)
        x = np.column_stack((xz_w, xm_w))
        A = np.zeros((nhost, 2, 2), dtype=np.float64)
        b = np.zeros((nhost, 2), dtype=np.float64)
        for i in range(nhost):
            mask = host_inverse == i
            A[i] = x[mask].T @ x[mask]
            b[i] = x[mask].T @ y_w[mask]
    else:
        bins = np.floor((xm - np.nanmin(xm)) / massbin_width).astype(np.int64)
        cells = np.rec.fromarrays((hosts, bins), names=("host", "massbin"))
        _, inverse, counts = np.unique(cells, return_inverse=True, return_counts=True)
        keep = counts[inverse] >= min_cell_n
        h = hosts[keep]
        z, yy, c = xz[keep], y[keep], inverse[keep]
        _, ci, cn = np.unique(c, return_inverse=True, return_counts=True)
        z_w = z - np.bincount(ci, weights=z)[ci] / cn[ci]
        y_w = yy - np.bincount(ci, weights=yy)[ci] / cn[ci]
        h_ids, hi = np.unique(h, return_inverse=True)
        A = np.bincount(hi, weights=z_w * z_w, minlength=len(h_ids))
        b = np.bincount(hi, weights=z_w * y_w, minlength=len(h_ids))
        nhost = len(h_ids)

    draws = rng.integers(0, nhost, size=(nboot, nhost))
    if massbin_width is None:
        boot_A = A[draws].sum(axis=1)
        boot_b = b[draws].sum(axis=1)
        coefs = np.full(nboot, math.nan, dtype=np.float64)
        for i in range(nboot):
            if np.linalg.matrix_rank(boot_A[i]) == 2:
                coefs[i] = np.linalg.solve(boot_A[i], boot_b[i])[0]
    else:
        boot_A = A[draws].sum(axis=1)
        boot_b = b[draws].sum(axis=1)
        coefs = boot_b / boot_A
        coefs = coefs[np.isfinite(coefs)]
    coefs = coefs[np.isfinite(coefs)]
    if not len(coefs):
        return math.nan, math.nan
    q = np.quantile(coefs, [0.025, 0.975])
    return float(q[0]), float(q[1])


def host_bootstrap_controls_interval(predictors, y, hosts, nboot, rng):
    """Bootstrap a host-fixed-effects coefficient with additional controls.

    Resampling whole hosts preserves the within-host dependence. The
    within-host cross-products are computed once, then summed with bootstrap
    multiplicities; repeated draws of one host represent independent bootstrap
    copies, each with its own fixed intercept.
    """
    X = np.column_stack(predictors).astype(np.float64, copy=False)
    y = np.asarray(y, dtype=np.float64)
    hosts = np.asarray(hosts)
    good = np.isfinite(X).all(axis=1) & np.isfinite(y)
    X, y, hosts = X[good], y[good], hosts[good]
    host_ids, host_inverse = np.unique(hosts, return_inverse=True)
    nhost, ncoef = len(host_ids), X.shape[1]
    if nhost < 2 or nboot <= 0 or len(y) <= ncoef:
        return math.nan, math.nan

    Xw = np.column_stack([within_host(X[:, j], hosts) for j in range(ncoef)])
    yw = within_host(y, hosts)
    cross = np.zeros((nhost, ncoef, ncoef), dtype=np.float64)
    response = np.zeros((nhost, ncoef), dtype=np.float64)
    for i in range(nhost):
        use = host_inverse == i
        cross[i] = Xw[use].T @ Xw[use]
        response[i] = Xw[use].T @ yw[use]

    draws = rng.integers(0, nhost, size=(nboot, nhost))
    boot_cross = cross[draws].sum(axis=1)
    boot_response = response[draws].sum(axis=1)
    coefficients = np.full(nboot, math.nan, dtype=np.float64)
    for i in range(nboot):
        if np.linalg.matrix_rank(boot_cross[i]) == ncoef:
            coefficients[i] = np.linalg.lstsq(
                boot_cross[i], boot_response[i], rcond=None
            )[0][0]
    coefficients = coefficients[np.isfinite(coefficients)]
    if not len(coefficients):
        return math.nan, math.nan
    q = np.quantile(coefficients, [0.025, 0.975])
    return float(q[0]), float(q[1])


def sigma3(v: np.ndarray) -> float:
    if len(v) < 2:
        return math.nan
    return float(np.sqrt(np.mean(np.var(v, axis=0, ddof=1))))


def host_fe_residual_scatter(predictors, y, hosts):
    """Residual scatter in log10(y) after host effects and listed controls."""
    X = np.column_stack(predictors).astype(np.float64, copy=False)
    y = np.asarray(y, dtype=np.float64)
    hosts = np.asarray(hosts)
    good = np.isfinite(X).all(axis=1) & np.isfinite(y)
    X, y, hosts = X[good], y[good], hosts[good]
    ngroup = len(np.unique(hosts))
    if len(y) <= ngroup + X.shape[1]:
        return math.nan
    Xw = np.column_stack([within_host(X[:, j], hosts) for j in range(X.shape[1])])
    yw = within_host(y, hosts)
    beta = np.linalg.lstsq(Xw, yw, rcond=None)[0]
    residual = yw - Xw @ beta
    dof = len(y) - ngroup - X.shape[1]
    return float(np.sqrt((residual @ residual) / dof))


def sigma_change_percent_per_metallicity_delta(slope, delta_logz=0.3):
    """Percent change in sigma for a fixed dex difference in metallicity."""
    if not np.isfinite(slope):
        return math.nan
    return float(100.0 * np.expm1(np.log(10.0) * slope * delta_logz))


def extract_rows(path: Path):
    rows = np.genfromtxt(path, delimiter=",", names=True, dtype=None, encoding="utf-8")
    return rows if rows.shape else rows.reshape(1)


def normalize_catalogue(path: Path):
    raw = extract_rows(path)
    names = set(raw.dtype.names or ())
    if "z_star" in names:
        return {
            "snapshot": raw["snapshot"],
            "redshift": raw["redshift"].astype(float),
            "z_star": raw["z_star"].astype(float),
            "z_star_number": np.full(len(raw), np.nan, dtype=float),
            "mstar_msun": raw["mstar_msun"].astype(float),
            "sfr100_msun_yr": raw["sfr100_msun_yr"].astype(float),
            "mgal_total_msun": np.full(len(raw), np.nan, dtype=float),
            "mdm_msun": (raw["mdm_msun"].astype(float) if "mdm_msun" in names
                         else np.full(len(raw), np.nan, dtype=float)),
            "rhalf3d_kpc": raw["rhalf3d_kpc"].astype(float),
            "age_mass_weighted_gyr": raw["age_mass_weighted_gyr"].astype(float),
            "vrot_global_kms": raw["vrot_global_kms"].astype(float),
            "host_id": raw["host_id"],
            "nstar": np.full(len(raw), -1, dtype=np.int64),
            "nstar_within_rhalf": np.full(len(raw), -1, dtype=np.int64),
            "particle_mstar_ratio": np.full(len(raw), np.nan, dtype=float),
            # Keep this legacy field as its own host-mass measure. Its contents
            # are not demonstrated to match the NewDD FoF DM+star mass, so do
            # not interpret the two catalogue adapters as mass-equivalent.
            "host_mass_msun": raw["host_mtot_msun"].astype(float),
            "vx_kms": raw["vx_kms"].astype(float),
            "vy_kms": raw["vy_kms"].astype(float),
            "vz_kms": raw["vz_kms"].astype(float),
            "sigma_global_kms": raw["sigma_global_kms"].astype(float),
            "sigma_mass_range_kms": raw["sigma_mass_range_kms"].astype(float),
            "sigma_within_rhalf_kms": raw["sigma_within_rhalf_kms"].astype(float),
            # The legacy ASCII table does not retain star-particle IDs, so
            # particle-level duplicate membership cannot be audited here.
            "duplicate_ok": np.ones(len(raw), dtype=bool),
            "duplicate_id_audit_available": False,
            "sigma_specs": [(key, key) for key in (
                "sigma_global_kms", "sigma_mass_range_kms", "sigma_within_rhalf_kms")],
            "source_host_mass_definition": "catalogue total host mass",
            "source_subhalo_mass_definition": "legacy catalogue Mdm field; definition unverified",
            "source_metallicity_measurement": (
                "legacy catalogue Z(star); weighting and absolute convention unverified"
            ),
        }

    required = {"is_central_max_stellar_mass",
                "mean_birth_metallicity_current_mass_weighted",
                "mean_birth_metallicity_number_weighted",
                "galaxy_stellar_mass_msun_h", "host_fof_dm_star_mass_msun_h",
                "galaxy_dm_star_mass_msun_h", "galaxy_dark_matter_mass_msun_h",
                "rhalf3d_stellar_mass_weighted_kpc",
                "sigma_star_1d_mass_weighted_within_rhalf_kms",
                "nstar_velocity_within_rhalf",
                "hubble_h", "host_index", "galaxy_vx_kms",
                "nstar_velocity_valid", "nstar_velocity_mass_valid"}
    required.add("stellar_particle_mass_sum_over_catalogue_mstar")
    if not required.issubset(names):
        raise ValueError(f"unsupported satellite catalogue fields: {sorted(names)}")
    raw = raw[raw["is_central_max_stellar_mass"] == 0]
    h = raw["hubble_h"].astype(float)
    duplicate_fields = {"nstar_id_duplicate", "nstar_ids_in_multi_galaxy_groups"}
    if not duplicate_fields.issubset(names):
        raise ValueError("NewGalFinder catalogue lacks the particle-ID duplicate audit fields")
    duplicate_ok = ((raw["nstar_id_duplicate"] == 0)
                    & (raw["nstar_ids_in_multi_galaxy_groups"] == 0))
    missing = np.full(len(raw), np.nan, dtype=float)
    return {
        "snapshot": raw["snapshot"],
        "redshift": raw["redshift"].astype(float),
        "z_star": raw["mean_birth_metallicity_current_mass_weighted"].astype(float),
        "z_star_number": raw["mean_birth_metallicity_number_weighted"].astype(float),
        "mstar_msun": raw["galaxy_stellar_mass_msun_h"].astype(float) / h,
        "sfr100_msun_yr": np.full(len(raw), np.nan, dtype=float),
        "mgal_total_msun": raw["galaxy_dm_star_mass_msun_h"].astype(float) / h,
        "mdm_msun": raw["galaxy_dark_matter_mass_msun_h"].astype(float) / h,
        "rhalf3d_kpc": raw["rhalf3d_stellar_mass_weighted_kpc"].astype(float),
        "age_mass_weighted_gyr": missing.copy(),
        "vrot_global_kms": missing.copy(),
        "host_id": raw["host_index"],
        "nstar": raw["nstar_velocity_mass_valid"].astype(np.int64),
        "nstar_within_rhalf": raw["nstar_velocity_within_rhalf"].astype(np.int64),
        "particle_mstar_ratio": raw[
            "stellar_particle_mass_sum_over_catalogue_mstar"].astype(float),
        "duplicate_ok": duplicate_ok,
        "duplicate_id_audit_available": True,
        "host_mass_msun": raw["host_fof_dm_star_mass_msun_h"].astype(float) / h,
        "vx_kms": raw["galaxy_vx_kms"].astype(float),
        "vy_kms": raw["galaxy_vy_kms"].astype(float),
        "vz_kms": raw["galaxy_vz_kms"].astype(float),
        "sigma_star_1d_mass_weighted_kms": raw["sigma_star_1d_mass_weighted_kms"].astype(float),
        "sigma_star_1d_mass_weighted_within_rhalf_kms": raw[
            "sigma_star_1d_mass_weighted_within_rhalf_kms"].astype(float),
        "sigma_star_1d_number_weighted_kms": raw["sigma_star_1d_number_weighted_kms"].astype(float),
        "sigma_specs": [
            ("sigma_star_1d_mass_weighted_kms", "mass-weighted 1D stellar dispersion"),
            ("sigma_star_1d_mass_weighted_within_rhalf_kms",
             "mass-weighted 1D stellar dispersion within 3D half-mass radius"),
            ("sigma_star_1d_number_weighted_kms", "number-weighted 1D stellar dispersion"),
        ],
        "source_host_mass_definition": (
            "FoF DM+star particle mass; gas and sinks omitted by NEWDD_SKIP_AMR"
        ),
        "source_subhalo_mass_definition": (
            "NewGalFinder-assigned DM+star particle mass; gas and sinks omitted"
        ),
        "source_metallicity_measurement": (
            "stellar-particle total birth-metallicity mass fraction zp; current-mass-weighted"
        ),
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("input", type=Path)
    ap.add_argument("output", type=Path)
    ap.add_argument("--min-satellite-mass", type=float, nargs="+", default=[1e8, 1e9])
    ap.add_argument("--stellar-mass-caliper-dex", type=float, default=0.2)
    ap.add_argument("--min-matched-pairs", type=int, default=4)
    ap.add_argument("--min-host-mass", type=float, default=0.0)
    ap.add_argument("--max-host-mass", type=float, default=float("inf"))
    ap.add_argument("--stellar-mass-bin-width-dex", type=float, default=0.3)
    ap.add_argument("--min-star-particles", type=int, default=100,
                    help="minimum stars with finite velocities and positive current masses")
    ap.add_argument("--min-aperture-star-particles", type=int, default=30,
                    help="minimum valid stellar velocities inside the 3D half-mass radius")
    ap.add_argument("--bootstrap", type=int, default=2000)
    args = ap.parse_args()
    if args.output.exists():
        raise FileExistsError(f"refusing to overwrite {args.output}")
    if args.min_aperture_star_particles < 2:
        raise ValueError("--min-aperture-star-particles must be at least 2")

    d = normalize_catalogue(args.input)
    z = d["z_star"]
    m = d["mstar_msun"]
    rhalf = d["rhalf3d_kpc"]
    age = d["age_mass_weighted_gyr"]
    vrot = d["vrot_global_kms"]
    host = d["host_id"]
    vel = np.column_stack((d["vx_kms"], d["vy_kms"], d["vz_kms"]))
    hm = d["host_mass_msun"]
    ssfr100 = d["sfr100_msun_yr"] / m
    sf_class_masks = {
        "star_forming": np.isfinite(ssfr100) & (ssfr100 >= 10.0**-10.5),
        "transitional": np.isfinite(ssfr100) & (ssfr100 > 0)
        & (ssfr100 > 10.0**-11.5) & (ssfr100 < 10.0**-10.5),
        "quiescent": np.isfinite(ssfr100) & (ssfr100 <= 10.0**-11.5),
    }
    nstar_ok = (d["nstar"] < 0) | (d["nstar"] >= args.min_star_particles)
    valid_internal_measurements = (
        np.isfinite(z) & (z > 0) & np.isfinite(m) & (m > 0)
    )
    valid_centre_velocity = np.isfinite(vel).all(axis=1)
    after_star_particle_cut = valid_internal_measurements & nstar_ok
    after_duplicate_cut = after_star_particle_cut & d["duplicate_ok"]
    after_host_mass_cut = (
        after_duplicate_cut & np.isfinite(hm)
        & (hm >= args.min_host_mass) & (hm < args.max_host_mass)
    )
    base = after_host_mass_cut

    summaries = []
    rng = np.random.default_rng(20261020)
    for mcut in args.min_satellite_mass:
        sel = base & (m >= mcut)
        particle_mass_ratio = d["particle_mstar_ratio"][sel]
        particle_mass_ratio = particle_mass_ratio[np.isfinite(particle_mass_ratio)]
        if particle_mass_ratio.size:
            particle_mass_abs_fractional_error = np.abs(particle_mass_ratio - 1.0)
            particle_mass_ratio_median = float(np.median(particle_mass_ratio))
            particle_mass_error_median_percent = float(
                100.0 * np.median(particle_mass_abs_fractional_error)
            )
            particle_mass_ratio_within_1pct = float(
                np.mean(particle_mass_abs_fractional_error <= 0.01)
            )
            particle_mass_ratio_within_5pct = float(
                np.mean(particle_mass_abs_fractional_error <= 0.05)
            )
        else:
            particle_mass_ratio_median = math.nan
            particle_mass_error_median_percent = math.nan
            particle_mass_ratio_within_1pct = math.nan
            particle_mass_ratio_within_5pct = math.nan
        host_ratios, low_medians, high_medians, matched_pair_counts = [], [], [], []
        n_hosts_orbital_with_pairs = 0
        n_hosts_orbital_invalid_sigma = 0
        orbital_sel = sel & valid_centre_velocity
        for h in np.unique(host[orbital_sel]):
            ix = np.flatnonzero(orbital_sel & (host == h))
            if len(ix) < 8:
                continue
            order = ix[np.argsort(z[ix])]
            half = len(order) // 2
            lo, hi = order[:half], order[-half:]
            if len(lo) < args.min_matched_pairs or len(hi) < args.min_matched_pairs:
                continue
            logm_lo, logm_hi = np.log10(m[lo]), np.log10(m[hi])
            ilo = np.argsort(logm_lo)
            ihi = np.argsort(logm_hi)
            i = j = 0
            pairs = []
            while i < len(ilo) and j < len(ihi):
                delta = float(logm_lo[ilo[i]] - logm_hi[ihi[j]])
                if abs(delta) <= args.stellar_mass_caliper_dex:
                    pairs.append((lo[ilo[i]], hi[ihi[j]]))
                    i += 1
                    j += 1
                elif delta < -args.stellar_mass_caliper_dex:
                    i += 1
                else:
                    j += 1
            if len(pairs) < args.min_matched_pairs:
                continue
            n_hosts_orbital_with_pairs += 1
            matched_lo = np.asarray([pair[0] for pair in pairs], dtype=int)
            matched_hi = np.asarray([pair[1] for pair in pairs], dtype=int)
            slo, shi = sigma3(vel[matched_lo]), sigma3(vel[matched_hi])
            if not (np.isfinite(slo) and np.isfinite(shi) and slo > 0 and shi > 0):
                n_hosts_orbital_invalid_sigma += 1
                continue
            host_ratios.append(math.log10(shi / slo))
            matched_pair_counts.append(len(pairs))
            low_medians.append(float(np.median(np.log10(m[matched_lo]))))
            high_medians.append(float(np.median(np.log10(m[matched_hi]))))

        ratios = np.asarray(host_ratios, dtype=float)
        if len(ratios):
            medratio = float(np.median(ratios))
            if args.bootstrap and len(ratios) > 1:
                boots = np.empty(args.bootstrap)
                for b in range(args.bootstrap):
                    boots[b] = np.median(ratios[rng.integers(0, len(ratios), len(ratios))])
                ci = np.quantile(boots, [0.025, 0.975])
            else:
                ci = [math.nan, math.nan]
            mass_offset = float(np.median(np.asarray(high_medians) - np.asarray(low_medians)))
        else:
            medratio = float("nan")
            ci = [math.nan, math.nan]
            mass_offset = float("nan")

        for sigma_field, sigma_definition in d["sigma_specs"]:
            sig = d[sigma_field].astype(float)
            # The aperture particle-count threshold applies only to the
            # half-mass-radius dispersion. Global dispersions do not depend on
            # how many stars happen to fall inside that aperture.
            if sigma_field == "sigma_star_1d_mass_weighted_within_rhalf_kms":
                aperture_count_ok = (
                    (d["nstar_within_rhalf"] < 0)
                    | (d["nstar_within_rhalf"] >= args.min_aperture_star_particles)
                )
            else:
                aperture_count_ok = np.ones(len(sig), dtype=bool)
            kin = sel & aperture_count_ok & np.isfinite(sig) & (sig > 0)
            beta, se, nfit, nhost_fit = fixed_host_ols(
                np.log10(z[kin]), np.log10(m[kin]), np.log10(sig[kin]), host[kin]
            ) if kin.sum() >= 10 else ((math.nan, math.nan), (math.nan, math.nan), int(kin.sum()), 0)
            z_number = d["z_star_number"]
            kin_number_z = kin & np.isfinite(z_number) & (z_number > 0)
            if kin_number_z.sum() >= 10:
                number_z_beta, number_z_se, number_z_n, number_z_nhost = fixed_host_ols(
                    np.log10(z_number[kin_number_z]), np.log10(m[kin_number_z]),
                    np.log10(sig[kin_number_z]), host[kin_number_z]
                )
                number_z_boot_lo, number_z_boot_hi = host_bootstrap_interval(
                    np.log10(z_number[kin_number_z]), np.log10(m[kin_number_z]),
                    np.log10(sig[kin_number_z]), host[kin_number_z],
                    args.bootstrap, rng,
                )
            else:
                number_z_beta, number_z_se = (math.nan, math.nan), (math.nan, math.nan)
                number_z_n, number_z_nhost = int(kin_number_z.sum()), 0
                number_z_boot_lo = number_z_boot_hi = math.nan
            sf_class_results = {}
            for sf_label, sf_mask in sf_class_masks.items():
                sf_kin = kin & sf_mask
                if sf_kin.sum() >= 10:
                    sf_beta, sf_se, sf_n, sf_nhost = fixed_host_ols(
                        np.log10(z[sf_kin]), np.log10(m[sf_kin]),
                        np.log10(sig[sf_kin]), host[sf_kin]
                    )
                    sf_boot_lo, sf_boot_hi = host_bootstrap_interval(
                        np.log10(z[sf_kin]), np.log10(m[sf_kin]),
                        np.log10(sig[sf_kin]), host[sf_kin], args.bootstrap, rng
                    )
                    sf_class_results.update({
                        f"internal_sigma_slope_{sf_label}": float(sf_beta[0]),
                        f"internal_sigma_slope_{sf_label}_se": float(sf_se[0]),
                        f"internal_sigma_slope_{sf_label}_bootstrap_95_low": sf_boot_lo,
                        f"internal_sigma_slope_{sf_label}_bootstrap_95_high": sf_boot_hi,
                        f"n_satellites_{sf_label}": sf_n,
                        f"n_hosts_{sf_label}": sf_nhost,
                    })
                else:
                    sf_class_results.update({
                        f"internal_sigma_slope_{sf_label}": math.nan,
                        f"internal_sigma_slope_{sf_label}_se": math.nan,
                        f"internal_sigma_slope_{sf_label}_bootstrap_95_low": math.nan,
                        f"internal_sigma_slope_{sf_label}_bootstrap_95_high": math.nan,
                        f"n_satellites_{sf_label}": int(sf_kin.sum()),
                        f"n_hosts_{sf_label}": 0,
                    })
            # Legacy catalogues also provide mean stellar age and half-mass
            # radius. These conditional fits are diagnostic: they show whether
            # the Z--sigma residual is shared with stellar-population age or
            # galaxy structure. NewGalFinder extracts lack these quantities.
            logz = np.log10(z[kin])
            logm = np.log10(m[kin])
            logsig = np.log10(sig[kin])
            mstar_only_residual_scatter = host_fe_residual_scatter(
                [logz, logm], logsig, host[kin]
            ) if kin.sum() >= 10 else math.nan
            kin_age = kin & np.isfinite(age) & (age > 0)
            sigma_age_beta = size_beta_control = joint_beta = np.full(2, math.nan)
            sigma_age_se = size_se_control = joint_se = np.full(2, math.nan)
            age_n = age_nhost = size_n = size_nhost = joint_n = joint_nhost = 0
            size_boot_lo = size_boot_hi = size_residual_scatter = math.nan
            if kin_age.sum() >= 10:
                b, e, age_n, age_nhost = fixed_host_ols_controls(
                    [np.log10(z[kin_age]), np.log10(m[kin_age]), np.log10(age[kin_age])],
                    np.log10(sig[kin_age]), host[kin_age]
                )
                sigma_age_beta, sigma_age_se = b[:2], e[:2]
            kin_size = kin & np.isfinite(rhalf) & (rhalf > 0)
            if kin_size.sum() >= 10:
                size_predictors = [
                    np.log10(z[kin_size]), np.log10(m[kin_size]),
                    np.log10(rhalf[kin_size]),
                ]
                b, e, size_n, size_nhost = fixed_host_ols_controls(
                    size_predictors,
                    np.log10(sig[kin_size]), host[kin_size]
                )
                size_beta_control, size_se_control = b[:2], e[:2]
                size_boot_lo, size_boot_hi = host_bootstrap_controls_interval(
                    size_predictors, np.log10(sig[kin_size]), host[kin_size],
                    args.bootstrap, rng,
                )
                size_residual_scatter = host_fe_residual_scatter(
                    size_predictors, np.log10(sig[kin_size]), host[kin_size],
                )
            kin_age_size = kin & np.isfinite(age) & (age > 0) & np.isfinite(rhalf) & (rhalf > 0)
            if kin_age_size.sum() >= 10:
                b, e, joint_n, joint_nhost = fixed_host_ols_controls(
                    [np.log10(z[kin_age_size]), np.log10(m[kin_age_size]),
                     np.log10(age[kin_age_size]), np.log10(rhalf[kin_age_size])],
                    np.log10(sig[kin_age_size]), host[kin_age_size]
                )
                joint_beta, joint_se = b[:2], e[:2]
            kin_total_mass = kin & np.isfinite(d["mgal_total_msun"]) & (d["mgal_total_msun"] > 0)
            total_mass_beta = total_mass_se = np.full(2, math.nan)
            total_mass_boot_lo = total_mass_boot_hi = math.nan
            total_mass_residual_scatter = math.nan
            total_mass_n = total_mass_nhost = 0
            if kin_total_mass.sum() >= 10:
                total_mass_predictors = [
                    np.log10(z[kin_total_mass]), np.log10(m[kin_total_mass]),
                    np.log10(d["mgal_total_msun"][kin_total_mass]),
                ]
                b, e, total_mass_n, total_mass_nhost = fixed_host_ols_controls(
                    total_mass_predictors,
                    np.log10(sig[kin_total_mass]), host[kin_total_mass]
                )
                total_mass_beta, total_mass_se = b[:2], e[:2]
                total_mass_residual_scatter = host_fe_residual_scatter(
                    total_mass_predictors, np.log10(sig[kin_total_mass]),
                    host[kin_total_mass],
                )
                total_mass_boot_lo, total_mass_boot_hi = host_bootstrap_controls_interval(
                    total_mass_predictors, np.log10(sig[kin_total_mass]),
                    host[kin_total_mass], args.bootstrap, rng,
                )
            kin_dark_matter = kin & np.isfinite(d["mdm_msun"]) & (d["mdm_msun"] > 0)
            dark_matter_beta = dark_matter_se = np.full(2, math.nan)
            dark_matter_boot_lo = dark_matter_boot_hi = math.nan
            dark_matter_residual_scatter = math.nan
            dark_matter_n = dark_matter_nhost = 0
            if kin_dark_matter.sum() >= 10:
                dark_matter_predictors = [
                    np.log10(z[kin_dark_matter]), np.log10(m[kin_dark_matter]),
                    np.log10(d["mdm_msun"][kin_dark_matter]),
                ]
                b, e, dark_matter_n, dark_matter_nhost = fixed_host_ols_controls(
                    dark_matter_predictors,
                    np.log10(sig[kin_dark_matter]), host[kin_dark_matter]
                )
                dark_matter_beta, dark_matter_se = b[:2], e[:2]
                dark_matter_residual_scatter = host_fe_residual_scatter(
                    dark_matter_predictors, np.log10(sig[kin_dark_matter]),
                    host[kin_dark_matter],
                )
                dark_matter_boot_lo, dark_matter_boot_hi = host_bootstrap_controls_interval(
                    dark_matter_predictors, np.log10(sig[kin_dark_matter]),
                    host[kin_dark_matter], args.bootstrap, rng,
                )
            boot_lo, boot_hi = host_bootstrap_interval(
                np.log10(z[kin]), np.log10(m[kin]), np.log10(sig[kin]), host[kin],
                args.bootstrap, rng,
            ) if kin.sum() >= 10 else (math.nan, math.nan)
            bbin, sebin, nbin, nhost_bin = fixed_host_massbin_slope(
                np.log10(z[kin]), np.log10(m[kin]), np.log10(sig[kin]), host[kin],
                width=args.stellar_mass_bin_width_dex,
            ) if kin.sum() >= 10 else (math.nan, math.nan, int(kin.sum()), 0)
            boot_bin_lo, boot_bin_hi = host_bootstrap_interval(
                np.log10(z[kin]), np.log10(m[kin]), np.log10(sig[kin]), host[kin],
                args.bootstrap, rng, massbin_width=args.stellar_mass_bin_width_dex,
            ) if kin.sum() >= 10 else (math.nan, math.nan)
            size_sel = sel & np.isfinite(rhalf) & (rhalf > 0)
            size_beta, size_se, _, _ = fixed_host_ols(
                np.log10(z[size_sel]), np.log10(m[size_sel]), np.log10(rhalf[size_sel]),
                host[size_sel],
            ) if size_sel.sum() >= 10 else ((math.nan, math.nan), (math.nan, math.nan), 0, 0)
            age_sel = sel & np.isfinite(age) & (age > 0)
            age_beta, age_se, _, _ = fixed_host_ols(
                np.log10(z[age_sel]), np.log10(m[age_sel]), np.log10(age[age_sel]),
                host[age_sel],
            ) if age_sel.sum() >= 10 else ((math.nan, math.nan), (math.nan, math.nan), 0, 0)
            rot_sel = kin & np.isfinite(vrot) & (vrot > 0)
            log_vrot_sigma = np.log10(vrot[rot_sel] / sig[rot_sel])
            rot_beta, rot_se, _, _ = fixed_host_ols(
                np.log10(z[rot_sel]), np.log10(m[rot_sel]), log_vrot_sigma, host[rot_sel]
            ) if rot_sel.sum() >= 10 else ((math.nan, math.nan), (math.nan, math.nan), 0, 0)
            summary_row = {
                "snapshot": int(d["snapshot"][0]),
                "redshift": float(d["redshift"][0]),
                "stellar_mass_cut_msun": mcut,
                "host_mass_definition": d["source_host_mass_definition"],
                "subhalo_mass_definition": d["source_subhalo_mass_definition"],
                "metallicity_measurement": d["source_metallicity_measurement"],
                "host_mass_min_msun": args.min_host_mass,
                "host_mass_max_msun": args.max_host_mass,
                "min_star_particles": args.min_star_particles,
                "internal_sigma_definition": sigma_definition,
                "n_satellites_in_catalogue": int(len(z)),
                "n_hosts_with_satellite_rows": int(len(np.unique(host))),
                "n_satellites_after_valid_internal_measurement_cut": int(
                    valid_internal_measurements.sum()
                ),
                "n_hosts_after_valid_internal_measurement_cut": int(
                    len(np.unique(host[valid_internal_measurements]))
                ),
                "n_satellites_after_star_particle_cut": int(after_star_particle_cut.sum()),
                "n_hosts_after_star_particle_cut": int(
                    len(np.unique(host[after_star_particle_cut]))
                ),
                "n_satellites_after_duplicate_id_cut": int(after_duplicate_cut.sum()),
                "n_hosts_after_duplicate_id_cut": int(
                    len(np.unique(host[after_duplicate_cut]))
                ),
                "n_satellites_after_host_mass_cut": int(after_host_mass_cut.sum()),
                "n_hosts_after_host_mass_cut": int(
                    len(np.unique(host[after_host_mass_cut]))
                ),
                "n_satellites_after_stellar_mass_cut": int(sel.sum()),
                "n_hosts_after_stellar_mass_cut": int(len(np.unique(host[sel]))),
                "n_satellites_selected": int(sel.sum()),
                "n_hosts_selected": int(len(np.unique(host[sel]))),
                "n_satellites_with_centre_velocity_for_orbital_check": int(
                    orbital_sel.sum()
                ),
                "n_satellites_with_particle_mass_consistency_check": int(
                    np.isfinite(d["particle_mstar_ratio"][sel]).sum()
                ),
                "median_particle_mass_sum_over_catalogue_mstar": particle_mass_ratio_median,
                "median_absolute_particle_mass_fractional_difference_percent":
                    particle_mass_error_median_percent,
                "fraction_particle_mass_sum_within_1pct_of_catalogue_mstar":
                    particle_mass_ratio_within_1pct,
                "fraction_particle_mass_sum_within_5pct_of_catalogue_mstar":
                    particle_mass_ratio_within_5pct,
                "n_satellites_excluded_duplicate_star_ids": int(
                    (nstar_ok & ~d["duplicate_ok"] & np.isfinite(z) & (z > 0)
                     & np.isfinite(m) & (m >= mcut) & np.isfinite(hm)
                     & (hm >= args.min_host_mass) & (hm < args.max_host_mass)).sum()
                ),
                "duplicate_star_id_audit_available": d["duplicate_id_audit_available"],
                "n_satellites_with_internal_sigma": int(kin.sum()),
                "n_satellites_excluded_low_aperture_particle_count": int(
                    (sel & ~aperture_count_ok).sum()
                ),
                "min_aperture_star_particles": args.min_aperture_star_particles,
                "n_hosts_internal_model": nhost_fit,
                "internal_sigma_log10Z_slope_at_fixed_logMstar_host": float(beta[0]),
                "cluster_robust_se": float(se[0]),
                "sigma_change_pct_per_0p3dex_Z_at_fixed_logMstar_host":
                    sigma_change_percent_per_metallicity_delta(beta[0]),
                "host_bootstrap_sigma_change_pct_low_per_0p3dex_Z":
                    sigma_change_percent_per_metallicity_delta(boot_lo),
                "host_bootstrap_sigma_change_pct_high_per_0p3dex_Z":
                    sigma_change_percent_per_metallicity_delta(boot_hi),
                "within_host_residual_scatter_dex_at_fixed_logMstar":
                    mstar_only_residual_scatter,
                "internal_sigma_log10Z_number_weighted_slope_at_fixed_logMstar_host": float(number_z_beta[0]),
                "cluster_robust_se_number_weighted_Z": float(number_z_se[0]),
                "sigma_change_pct_per_0p3dex_number_weighted_Z":
                    sigma_change_percent_per_metallicity_delta(number_z_beta[0]),
                "host_bootstrap_95_low_number_weighted_Z": number_z_boot_lo,
                "host_bootstrap_95_high_number_weighted_Z": number_z_boot_hi,
                "host_bootstrap_sigma_change_pct_low_number_weighted_Z":
                    sigma_change_percent_per_metallicity_delta(number_z_boot_lo),
                "host_bootstrap_sigma_change_pct_high_number_weighted_Z":
                    sigma_change_percent_per_metallicity_delta(number_z_boot_hi),
                "n_satellites_number_weighted_Z_model": number_z_n,
                "n_hosts_number_weighted_Z_model": number_z_nhost,
                "internal_sigma_slope_controlling_log_age": float(sigma_age_beta[0]),
                "internal_sigma_slope_controlling_log_age_se": float(sigma_age_se[0]),
                "n_satellites_age_control_model": age_n,
                "internal_sigma_slope_controlling_log_rhalf": float(size_beta_control[0]),
                "internal_sigma_slope_controlling_log_rhalf_se": float(size_se_control[0]),
                "sigma_change_pct_per_0p3dex_Z_controlling_log_rhalf":
                    sigma_change_percent_per_metallicity_delta(size_beta_control[0]),
                "host_bootstrap_95_low_controlling_log_rhalf": size_boot_lo,
                "host_bootstrap_95_high_controlling_log_rhalf": size_boot_hi,
                "host_bootstrap_sigma_change_pct_low_controlling_log_rhalf":
                    sigma_change_percent_per_metallicity_delta(size_boot_lo),
                "host_bootstrap_sigma_change_pct_high_controlling_log_rhalf":
                    sigma_change_percent_per_metallicity_delta(size_boot_hi),
                "within_host_residual_scatter_dex_controlling_log_rhalf":
                    size_residual_scatter,
                "n_satellites_rhalf_control_model": size_n,
                "internal_sigma_slope_controlling_log_age_and_log_rhalf": float(joint_beta[0]),
                "internal_sigma_slope_controlling_log_age_and_log_rhalf_se": float(joint_se[0]),
                "n_satellites_age_rhalf_control_model": joint_n,
                "internal_sigma_slope_controlling_log_subhalo_dm_star_mass": float(total_mass_beta[0]),
                "internal_sigma_slope_controlling_log_subhalo_dm_star_mass_se": float(total_mass_se[0]),
                "sigma_change_pct_per_0p3dex_Z_controlling_subhalo_dm_star_mass":
                    sigma_change_percent_per_metallicity_delta(total_mass_beta[0]),
                "host_bootstrap_sigma_change_pct_low_controlling_subhalo_dm_star_mass":
                    sigma_change_percent_per_metallicity_delta(total_mass_boot_lo),
                "host_bootstrap_sigma_change_pct_high_controlling_subhalo_dm_star_mass":
                    sigma_change_percent_per_metallicity_delta(total_mass_boot_hi),
                "within_host_residual_scatter_dex_controlling_subhalo_dm_star_mass":
                    total_mass_residual_scatter,
                "host_bootstrap_95_low_subhalo_dm_star_mass_control": total_mass_boot_lo,
                "host_bootstrap_95_high_subhalo_dm_star_mass_control": total_mass_boot_hi,
                "n_satellites_subhalo_dm_star_mass_control_model": total_mass_n,
                "n_hosts_subhalo_dm_star_mass_control_model": total_mass_nhost,
                "internal_sigma_slope_controlling_log_subhalo_dark_matter_mass": float(dark_matter_beta[0]),
                "internal_sigma_slope_controlling_log_subhalo_dark_matter_mass_se": float(dark_matter_se[0]),
                "sigma_change_pct_per_0p3dex_Z_controlling_subhalo_dark_matter_mass":
                    sigma_change_percent_per_metallicity_delta(dark_matter_beta[0]),
                "host_bootstrap_sigma_change_pct_low_controlling_subhalo_dark_matter_mass":
                    sigma_change_percent_per_metallicity_delta(dark_matter_boot_lo),
                "host_bootstrap_sigma_change_pct_high_controlling_subhalo_dark_matter_mass":
                    sigma_change_percent_per_metallicity_delta(dark_matter_boot_hi),
                "within_host_residual_scatter_dex_controlling_subhalo_dark_matter_mass":
                    dark_matter_residual_scatter,
                "host_bootstrap_95_low_subhalo_dark_matter_mass_control": dark_matter_boot_lo,
                "host_bootstrap_95_high_subhalo_dark_matter_mass_control": dark_matter_boot_hi,
                "n_satellites_subhalo_dark_matter_mass_control_model": dark_matter_n,
                "n_hosts_subhalo_dark_matter_mass_control_model": dark_matter_nhost,
                "host_bootstrap_95_low": boot_lo,
                "host_bootstrap_95_high": boot_hi,
                "slope_at_fixed_host_and_stellar_mass_bins": bbin,
                "cluster_robust_se_mass_bins": sebin,
                "host_bootstrap_95_low_mass_bins": boot_bin_lo,
                "host_bootstrap_95_high_mass_bins": boot_bin_hi,
                "mass_bin_width_dex": args.stellar_mass_bin_width_dex,
                "n_satellites_mass_bin_model": nbin,
                "n_hosts_mass_bin_model": nhost_bin,
                "rhalf3d_logZ_slope_at_fixed_logMstar_host": float(size_beta[0]),
                "rhalf3d_cluster_robust_se": float(size_se[0]),
                "age_logZ_slope_at_fixed_logMstar_host": float(age_beta[0]),
                "age_cluster_robust_se": float(age_se[0]),
                "vrot_over_sigma_logZ_slope_at_fixed_logMstar_host": float(rot_beta[0]),
                "vrot_over_sigma_cluster_robust_se": float(rot_se[0]),
                "spearman_logZ_log_sigma_raw": corr(np.log10(z[kin]), np.log10(sig[kin])),
                "n_hosts_orbital_with_mass_matched_pairs": n_hosts_orbital_with_pairs,
                "n_hosts_orbital_excluded_nonpositive_sigma": n_hosts_orbital_invalid_sigma,
                "n_hosts_orbital_split": len(ratios),
                "n_matched_satellite_pairs_orbital_total": int(sum(matched_pair_counts)),
                "median_matched_pairs_per_orbital_host": (
                    float(np.median(matched_pair_counts)) if matched_pair_counts else math.nan
                ),
                "minimum_matched_pairs_per_orbital_host": (
                    int(min(matched_pair_counts)) if matched_pair_counts else 0
                ),
                "mass_match_caliper_dex": args.stellar_mass_caliper_dex,
                "min_matched_pairs_per_host": args.min_matched_pairs,
                "median_log10_orbital_sigma_highZ_over_lowZ": medratio,
                "bootstrap_95_low": float(ci[0]),
                "bootstrap_95_high": float(ci[1]),
                "median_log10Mstar_highZ_minus_lowZ": mass_offset,
            }
            summary_row.update(sf_class_results)
            summaries.append(summary_row)

    temp_output = args.output.with_name(
        f".{args.output.name}.partial.{os.getpid()}"
    )
    if temp_output.exists():
        raise FileExistsError(f"refusing to overwrite {temp_output}")
    with temp_output.open("x", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(summaries[0]))
        writer.writeheader()
        writer.writerows(summaries)
    os.link(temp_output, args.output)
    temp_output.unlink()
    for row in summaries:
        print(row)
    print(f"summary={args.output}")


if __name__ == "__main__":
    main()
