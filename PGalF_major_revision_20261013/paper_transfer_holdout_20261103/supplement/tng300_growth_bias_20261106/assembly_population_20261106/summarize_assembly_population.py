#!/usr/bin/env python3
"""Create branch-paired summaries for the TNG300 assembly/population analysis.

This script reads the per-branch and per-view products of
``analyze_assembly_population.py``. It does not read simulation cutouts or
recompute particle-level measurements. Projections are paired before their
branch median is taken; bootstrap resampling then uses the 36 branches.
"""
from __future__ import annotations

import hashlib
import argparse
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
INPUT = HERE / "results"
OUTPUT = HERE / "results"
BOOTSTRAPS = 20_000
SEED = 20261106
METHOD_LABELS = {
    "isophotal_particle_selection": "isophotal",
    "b010_linked_group": "FoF b=0.10",
    "b012_linked_group": "FoF b=0.12",
    "b020_linked_group": "FoF b=0.20",
}
RADIUS_COLUMNS = {
    "mass": "R50_current_stellar_mass_pkpc",
    "full": "R50_full_stellar_ages_and_metallicities_pkpc",
    "uniform": "R50_uniform_age_and_metallicity_pkpc",
    "uniform_age_varying_Z": "R50_uniform_age_varying_metallicity_pkpc",
    "varying_age_uniform_Z": "R50_varying_age_uniform_metallicity_pkpc",
}


def stable_seed(key: str) -> int:
    digest = hashlib.sha256(key.encode("utf-8")).digest()
    return SEED + int.from_bytes(digest[:4], "little")


def summarize(values, *, key: str, unit: str, metric: str, band: str = "",
              condition: str = "", interval: str = "", notes: str = ""):
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if not values.size:
        raise ValueError(f"No finite values for {key}")
    rng = np.random.default_rng(stable_seed(key))
    draws = rng.integers(0, values.size, size=(BOOTSTRAPS, values.size))
    boot_medians = np.median(values[draws], axis=1)
    lo, hi = np.percentile(boot_medians, [2.5, 97.5])
    p16, p84 = np.percentile(values, [16, 84])
    return {
        "metric": metric,
        "band": band,
        "condition": condition,
        "interval": interval,
        "unit": unit,
        "n_branches": int(values.size),
        "median": float(np.median(values)),
        "p16": float(p16),
        "p84": float(p84),
        "median_ci95_low": float(lo),
        "median_ci95_high": float(hi),
        "notes": notes,
    }


def branch_medians(frame: pd.DataFrame, value: str, groups: list[str]):
    return frame.groupby(groups, as_index=False)[value].median()


def main():
    global INPUT, OUTPUT
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, default=HERE / "results",
                        help="directory containing the analysis CSV products")
    parser.add_argument("--output", type=Path, default=None,
                        help="directory for compact summaries (defaults to --input)")
    args = parser.parse_args()
    INPUT = args.input.resolve()
    OUTPUT = (args.output if args.output is not None else INPUT).resolve()
    OUTPUT.mkdir(parents=True, exist_ok=True)

    rows = []
    budget = pd.read_csv(INPUT / "assembly_mass_budget_per_branch.csv")
    long = budget[(budget.snapshot_start == 50) & (budget.snapshot_end == 84)].copy()
    budget_metrics = {
        "net central-subhalo mass change": ("delta_over_start", 1.0),
        "persistent-star current-mass change": ("persistent_change_over_start", 1.0),
        "pre-existing stars newly assigned to central": ("acquired_old_over_start", 1.0),
        "stars formed after snapshot 50, central at snapshot 84": ("new_formed_over_start", 1.0),
        "snapshot-50 central stars absent from snapshot-84 central": ("lost_over_start", -1.0),
    }
    for metric, (column, sign) in budget_metrics.items():
        rows.append(summarize(
            sign * long[column], key=f"budget:{metric}", unit="M_start",
            metric=metric, interval="snapshot 50 to 84 (z=0.997 to 0.197)",
            notes="Each branch is normalized by its snapshot-50 central stellar mass."))
    sat_fraction = (
        long.mass_acquired_old_start_satellite_mask_msun
        / long.acquired_preexisting_mass_msun
    )
    rows.append(summarize(
        sat_fraction, key="budget:resolved-satellite-source-fraction",
        unit="fraction of pre-existing newly central mass",
        metric="Earlier resolved-satellite-mask fraction of pre-existing newly central mass",
        interval="snapshot 50 to 84 (z=0.997 to 0.197)",
        notes="Particle mass in the earlier RemoveSatelliteGalaxy mask divided by all pre-existing mass newly assigned to the central."))

    # Branch-paired capture fractions at the final output.
    capture = pd.read_csv(INPUT / "selection_population_capture_per_view.csv")
    capture = capture[
        (capture.snapshot == 84)
        & capture.method.isin(METHOD_LABELS)
    ].copy()
    branch_capture = branch_medians(
        capture, "capture_fraction_of_class",
        ["sample_rank", "band", "method", "population_class"])
    for keys, group in branch_capture.groupby(["band", "method", "population_class"]):
        band, method, pop_class = keys
        rows.append(summarize(
            group.capture_fraction_of_class,
            key=f"capture:{band}:{method}:{pop_class}", unit="fraction",
            metric=f"capture: {pop_class}; {METHOD_LABELS[method]}", band=band,
            interval="snapshot 84 (z=0.197)",
            notes="Median over three projections per branch, followed by a branch median."))

    # Size growth from paired endpoint measurements, retaining projection pairing.
    sizes = pd.read_csv(INPUT / "central_size_metrics_per_view.csv")
    sizes = sizes[sizes.snapshot.isin([50, 84])].copy()
    key_cols = ["sample_rank", "band", "observer_condition", "line_of_sight"]
    wide = sizes.pivot(index=key_cols, columns="snapshot", values=list(RADIUS_COLUMNS.values()))
    growth_rows = []
    for short, col in RADIUS_COLUMNS.items():
        growth = np.log10(wide[(col, 84)] / wide[(col, 50)])
        growth_rows.append(pd.DataFrame({
            **{name: wide.index.get_level_values(name) for name in key_cols},
            "quantity": short,
            "growth_dex": growth.to_numpy(),
        }))
    growth = pd.concat(growth_rows, ignore_index=True)
    paired = growth.pivot(
        index=key_cols, columns="quantity", values="growth_dex").reset_index()
    paired["combined_population_gradient_effect"] = paired.full - paired.uniform
    paired["age_shapley_effect"] = 0.5 * (
        (paired.varying_age_uniform_Z - paired.uniform)
        + (paired.full - paired.uniform_age_varying_Z)
    )
    paired["metallicity_shapley_effect"] = 0.5 * (
        (paired.uniform_age_varying_Z - paired.uniform)
        + (paired.full - paired.varying_age_uniform_Z)
    )
    branch_size = paired.groupby(
        ["sample_rank", "band", "observer_condition"], as_index=False
    ).median(numeric_only=True)
    growth_names = {
        "mass": "central half-mass-radius growth",
        "full": "full age-and-metallicity-weighted half-light-radius growth",
        "uniform": "uniform-population half-light-radius growth",
        "combined_population_gradient_effect": "paired age-and-metallicity gradient effect on size growth",
        "age_shapley_effect": "age contribution (Shapley split)",
        "metallicity_shapley_effect": "metallicity contribution (Shapley split)",
    }
    for (band, condition), group in branch_size.groupby(["band", "observer_condition"]):
        for col, metric in growth_names.items():
            rows.append(summarize(
                group[col], key=f"size:{band}:{condition}:{col}", unit="dex",
                metric=metric, band=band, condition=condition,
                interval="snapshot 50 to 84 (z=0.997 to 0.197)",
                notes="Paired within branch and projection; then median over projections and branches. Age and metallicity Shapley medians need not sum to the combined median."))
    paired_conditions = paired.pivot(
        index=["sample_rank", "band", "line_of_sight"],
        columns="observer_condition", values="full").reset_index()
    paired_conditions["physical_minus_fixed"] = (
        paired_conditions.physical_observer_redshift
        - paired_conditions.fixed_observer_z0p197
    )
    observer_branch = paired_conditions.groupby(
        ["sample_rank", "band"], as_index=False).physical_minus_fixed.median()
    for band, group in observer_branch.groupby("band"):
        rows.append(summarize(
            group.physical_minus_fixed, key=f"observer:physical-minus-fixed:{band}",
            unit="dex", metric="physical-observer minus fixed-observer light-size growth",
            band=band, interval="snapshot 50 to 84",
            notes="Paired by branch and projection; fixed observer redshift is 0.197284."))

    pd.DataFrame(rows).to_csv(OUTPUT / "key_results.csv", index=False, float_format="%.8g")
    paired.to_csv(OUTPUT / "central_size_growth_paired_per_view.csv", index=False,
                  float_format="%.10g")

    # Combine radial bins by stellar mass for <1 and >1 current-mass R50.
    profiles = pd.read_csv(INPUT / "radial_population_profiles_per_view.csv")
    profiles = profiles[profiles.observer_condition.isin([
        "physical_observer_redshift", "fixed_observer_z0p197"]
    )].copy()
    profiles["inner"] = profiles.radial_bin_R50mass.isin(["0-0.5", "0.5-1"])
    profiles["light_weight"] = profiles.current_stellar_mass_msun / (
        10.0 ** profiles.log10_current_mass_to_light_ratio
    )
    radial_rows = []
    radial_keys = ["sample_rank", "snapshot", "observer_condition", "band",
                   "line_of_sight", "inner"]
    for keys, group in profiles.groupby(radial_keys):
        rank, snap, condition, band, los, inner = keys
        mass = group.current_stellar_mass_msun.to_numpy()
        total_mass = mass.sum()
        if total_mass <= 0:
            continue
        radial_rows.append({
            "sample_rank": rank,
            "snapshot": snap,
            "observer_condition": condition,
            "band": band,
            "line_of_sight": los,
            "inner": inner,
            "mass_weighted_mean_age_gyr": np.average(
                group.mass_weighted_mean_age_gyr, weights=mass),
            "mass_weighted_mean_log10_metal_mass_fraction": np.average(
                group.mass_weighted_mean_log10_metal_mass_fraction, weights=mass),
            "log10_current_mass_to_light_ratio": np.log10(
                total_mass / group.light_weight.sum()),
        })
    radial = pd.DataFrame(radial_rows)
    radial_value_cols = [
        "mass_weighted_mean_age_gyr",
        "mass_weighted_mean_log10_metal_mass_fraction",
        "log10_current_mass_to_light_ratio",
    ]
    radial_wide = radial.pivot(
        index=["sample_rank", "snapshot", "observer_condition", "band", "line_of_sight"],
        columns="inner", values=radial_value_cols)
    radial_diffs = []
    for col in radial_value_cols:
        frame = radial_wide[(col, False)] - radial_wide[(col, True)]
        radial_diffs.append(frame.rename(col + "_outer_minus_inner"))
    radial_diff = pd.concat(radial_diffs, axis=1).reset_index()
    radial_branch = radial_diff.groupby(
        ["sample_rank", "snapshot", "observer_condition", "band"],
        as_index=False).median(numeric_only=True)
    radial_summary = []
    for (snap, condition, band), group in radial_branch.groupby(
        ["snapshot", "observer_condition", "band"]):
        for col in [name + "_outer_minus_inner" for name in radial_value_cols]:
            radial_summary.append(summarize(
                group[col], key=f"radial:{snap}:{condition}:{band}:{col}",
                unit="Gyr" if col.startswith("mass_weighted_mean_age") else "dex",
                metric=col.replace("_outer_minus_inner", " (R>R50mass minus R<R50mass)"),
                band=band, condition=condition,
                interval=f"snapshot {snap}",
                notes="Radial bins below/above one projected current-mass half-mass radius are combined within each view; views are reduced to branch medians."))
    pd.DataFrame(radial_summary).to_csv(
        OUTPUT / "radial_gradient_summary.csv", index=False, float_format="%.8g")
    radial_diff.to_csv(OUTPUT / "radial_gradients_per_view.csv", index=False,
                       float_format="%.10g")

    print(f"Wrote {OUTPUT / 'key_results.csv'}")
    print(f"Wrote {OUTPUT / 'radial_gradient_summary.csv'}")
    print(f"Wrote {OUTPUT / 'central_size_growth_paired_per_view.csv'}")
    print(f"Wrote {OUTPUT / 'radial_gradients_per_view.csv'}")


if __name__ == "__main__":
    main()
