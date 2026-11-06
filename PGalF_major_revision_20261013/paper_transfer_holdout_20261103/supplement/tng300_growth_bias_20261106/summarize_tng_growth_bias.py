#!/usr/bin/env python3
"""Summarize selection bias in TNG300-1 stellar mass and size growth."""
from __future__ import annotations

import csv
import json
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np

HERE = Path(__file__).resolve().parent
PAPER = HERE.parents[1]
SCORE_DIR = HERE / "growth_scores"
FIGURE = PAPER / "figures/tng300_tree_growth_bias.pdf"
SEED = 20261106
N_BOOT = 20000
RANKS = sorted(int(row["sample_rank"]) for row in csv.DictReader(
    (HERE / "lineage_cohort.csv").open()) if int(row["central_at_all_targets"]) == 1)
SNAPS = (84, 67, 50)
METHODS = (
    "isophotal_particle_selection",
    "b010_linked_group",
    "b020_linked_group",
    "projected_aperture_30_pkpc",
    "Darwin_calibrated_sphere",
)
INTERVALS = ((50, 67, "z1_to_z0p5"), (67, 84, "z0p5_to_z0p2"),
             (50, 84, "z1_to_z0p2"))
CONDITIONS = ("actual_observer_redshift", "fixed_observer_z0p197")
OBSERVER_REDSHIFTS = (0.19728418237601, 0.503047523244883, 0.99729422578194)
SNAP_TO_REDSHIFT = {84: 0.19728418237601, 67: 0.503047523244883, 50: 0.99729422578194}


def read_csv(path: Path):
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream))


def quantile_interval(values):
    q = np.quantile(np.asarray(values, dtype=float), [0.025, 0.975])
    return [float(q[0]), float(q[1])]


def bootstrap(values, rng: np.random.Generator, require_full_cohort: bool = True):
    values = np.asarray(values, dtype=float)
    if values.size < 1:
        raise RuntimeError("Cannot bootstrap an empty tree sample")
    if require_full_cohort and values.size != len(RANKS):
        raise RuntimeError(f"Expected {len(RANKS)} independent tree values, received {values.size}")
    sample_idx = rng.integers(0, values.size, size=(N_BOOT, values.size))
    draws = values[sample_idx]
    mean_ci = quantile_interval(np.mean(draws, axis=1))
    median_ci = quantile_interval(np.median(draws, axis=1))
    abs_draws = np.abs(draws)
    return {
        "n_trees": int(values.size),
        "mean": float(np.mean(values)),
        "mean_ci95": mean_ci,
        "median": float(np.median(values)),
        "median_ci95": median_ci,
        "mean_absolute_error": float(np.mean(np.abs(values))),
        "mean_absolute_error_ci95": quantile_interval(np.mean(abs_draws, axis=1)),
        "median_absolute_error": float(np.median(np.abs(values))),
        "median_absolute_error_ci95": quantile_interval(np.median(abs_draws, axis=1)),
    }


def write_csv(path: Path, rows):
    if path.exists():
        raise FileExistsError(f"Refusing to overwrite existing analysis output: {path.name}")
    if not rows:
        raise RuntimeError(f"No rows to write for {path.name}")
    with path.open("x", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main():
    if FIGURE.exists():
        raise FileExistsError(f"Refusing to overwrite existing figure: {FIGURE.name}")
    score_files = [SCORE_DIR / f"tree_rank_{rank:02d}_growth_scores.csv" for rank in RANKS]
    if any(not path.is_file() for path in score_files):
        missing = [p.name for p in score_files if not p.is_file()]
        raise FileNotFoundError(f"Missing per-tree score files: {missing[:5]}")

    # Keep the observer-redshift cases separate instead of pooling views or outputs.
    indexed = defaultdict(list)
    central_mass_audit = {}
    for path in score_files:
        rows = read_csv(path)
        if len(rows) != 324:
            raise RuntimeError(f"Expected 324 rows in {path.name}, found {len(rows)}")
        rank = int(path.stem.split("_")[2])
        for row in rows:
            key = (rank, int(row["snapshot"]), row["band"], row["method"], row["line_of_sight"])
            indexed[key].append(row)
            central_mass_audit[(rank, int(row["snapshot"]))] = float(row["central_catalog_particle_mass_ratio"])
    max_catalog_mass_log_error = max(abs(np.log10(v)) for v in central_mass_audit.values())

    def select_row(rank, snap, band, method, los, condition):
        candidates = indexed[(rank, snap, band, method, los)]
        if condition == "fixed_observer_z0p197":
            target_z = OBSERVER_REDSHIFTS[0]
        else:
            target_z = SNAP_TO_REDSHIFT[snap]
        matches = [r for r in candidates if abs(float(r["observer_redshift"]) - target_z) < 1e-8]
        if len(matches) != 1:
            raise RuntimeError(
                f"Expected one row for rank={rank}, snap={snap}, band={band}, "
                f"method={method}, view={los}, condition={condition}; got {len(matches)}"
            )
        return matches[0]

    per_tree = []
    for band in ("V", "g"):
        for condition in CONDITIONS:
            for method in METHODS:
                for rank in RANKS:
                    per_view = []
                    for los in ("x", "y", "z"):
                        rows = {snap: select_row(rank, snap, band, method, los, condition)
                                for snap in SNAPS}
                        refs = {snap: select_row(rank, snap, band,
                                                 "central_bound_stars_reference", los, condition)
                                for snap in SNAPS}
                        measurements = {}
                        for snap in SNAPS:
                            m = float(rows[snap]["stellar_mass_msun"])
                            mb = float(refs[snap]["central_bound_stellar_mass_msun"])
                            r50 = float(rows[snap]["projected_R50_pkpc"])
                            rb = float(refs[snap]["central_bound_R50_pkpc"])
                            if min(m, mb, r50, rb) <= 0:
                                raise RuntimeError(f"Non-positive mass or radius for rank={rank}, {method}")
                            measurements[snap] = (m, mb, r50, rb)
                        interval_biases = {}
                        for start, end, label in INTERVALS:
                            m_start, mb_start, r_start, rb_start = measurements[start]
                            m_end, mb_end, r_end, rb_end = measurements[end]
                            mass_growth = float(np.log10(m_end / m_start))
                            bound_mass_growth = float(np.log10(mb_end / mb_start))
                            radius_growth = float(np.log10(r_end / r_start))
                            bound_radius_growth = float(np.log10(rb_end / rb_start))
                            interval_biases[label] = {
                                "stellar_mass_growth_dex": mass_growth,
                                "bound_stellar_mass_growth_dex": bound_mass_growth,
                                "mass_growth_bias_dex": mass_growth - bound_mass_growth,
                                "R50_growth_dex": radius_growth,
                                "bound_R50_growth_dex": bound_radius_growth,
                                "R50_growth_bias_dex": radius_growth - bound_radius_growth,
                            }
                        mass_growth = interval_biases["z1_to_z0p2"]["stellar_mass_growth_dex"]
                        bound_mass_growth = interval_biases["z1_to_z0p2"]["bound_stellar_mass_growth_dex"]
                        radius_growth = interval_biases["z1_to_z0p2"]["R50_growth_dex"]
                        bound_radius_growth = interval_biases["z1_to_z0p2"]["bound_R50_growth_dex"]
                        per_view.append({
                            "line_of_sight": los,
                            **{f"{key}_{interval}": value
                               for interval, values in interval_biases.items()
                               for key, value in values.items()},
                            "mass_ratio_z0p2": measurements[84][0] / measurements[84][1],
                            "mass_ratio_z0p5": measurements[67][0] / measurements[67][1],
                            "mass_ratio_z1": measurements[50][0] / measurements[50][1],
                        })
                    view_medians = {}
                    for _start, _end, interval in INTERVALS:
                        for metric in ("stellar_mass_growth_dex", "bound_stellar_mass_growth_dex",
                                       "mass_growth_bias_dex", "R50_growth_dex",
                                       "bound_R50_growth_dex", "R50_growth_bias_dex"):
                            key = f"{metric}_{interval}"
                            view_medians[key] = float(np.median([v[key] for v in per_view]))
                    per_tree.append({
                        "sample_rank": rank,
                        "band": band,
                        "condition": condition,
                        "method": method,
                        **view_medians,
                        "mass_growth_bias_dex": view_medians["mass_growth_bias_dex_z1_to_z0p2"],
                        "R50_growth_bias_dex": view_medians["R50_growth_bias_dex_z1_to_z0p2"],
                        "stellar_mass_growth_dex": view_medians["stellar_mass_growth_dex_z1_to_z0p2"],
                        "bound_stellar_mass_growth_dex": view_medians[
                            "bound_stellar_mass_growth_dex_z1_to_z0p2"],
                        "R50_growth_dex": view_medians["R50_growth_dex_z1_to_z0p2"],
                        "bound_R50_growth_dex": view_medians["bound_R50_growth_dex_z1_to_z0p2"],
                        "median_mass_ratio_z0p2": float(np.median([v["mass_ratio_z0p2"] for v in per_view])),
                        "median_mass_ratio_z0p5": float(np.median([v["mass_ratio_z0p5"] for v in per_view])),
                        "median_mass_ratio_z1": float(np.median([v["mass_ratio_z1"] for v in per_view])),
                    })

    rng = np.random.default_rng(SEED)
    summary = []
    interval_summary = []
    summaries_by_key = {}
    for band in ("V", "g"):
        for condition in CONDITIONS:
            for method in METHODS:
                rows = [r for r in per_tree if r["band"] == band and r["condition"] == condition
                        and r["method"] == method]
                mass_summary = bootstrap(np.array([r["mass_growth_bias_dex"] for r in rows]), rng)
                size_summary = bootstrap(np.array([r["R50_growth_bias_dex"] for r in rows]), rng)
                summaries_by_key[(band, condition, method)] = (mass_summary, size_summary)
                summary.append({
                    "band": band, "condition": condition, "method": method,
                    "n_descendant_trees": len(rows),
                    "mean_mass_growth_bias_dex": mass_summary["mean"],
                    "mean_mass_growth_bias_ci95": json.dumps(mass_summary["mean_ci95"]),
                    "median_mass_growth_bias_dex": mass_summary["median"],
                    "median_mass_growth_bias_ci95": json.dumps(mass_summary["median_ci95"]),
                    "mean_absolute_mass_growth_bias_dex": mass_summary["mean_absolute_error"],
                    "mean_absolute_mass_growth_bias_ci95": json.dumps(mass_summary["mean_absolute_error_ci95"]),
                    "mean_R50_growth_bias_dex": size_summary["mean"],
                    "mean_R50_growth_bias_ci95": json.dumps(size_summary["mean_ci95"]),
                    "median_R50_growth_bias_dex": size_summary["median"],
                    "median_R50_growth_bias_ci95": json.dumps(size_summary["median_ci95"]),
                    "mean_absolute_R50_growth_bias_dex": size_summary["mean_absolute_error"],
                    "mean_absolute_R50_growth_bias_ci95": json.dumps(size_summary["mean_absolute_error_ci95"]),
                })
                for _start, _end, interval in INTERVALS:
                    interval_mass = bootstrap(np.array([
                        r[f"mass_growth_bias_dex_{interval}"] for r in rows]), rng)
                    interval_size = bootstrap(np.array([
                        r[f"R50_growth_bias_dex_{interval}"] for r in rows]), rng)
                    interval_summary.append({
                        "band": band,
                        "condition": condition,
                        "method": method,
                        "interval": interval,
                        "n_descendant_trees": len(rows),
                        "mean_mass_growth_bias_dex": interval_mass["mean"],
                        "mean_mass_growth_bias_ci95": json.dumps(interval_mass["mean_ci95"]),
                        "median_mass_growth_bias_dex": interval_mass["median"],
                        "median_mass_growth_bias_ci95": json.dumps(interval_mass["median_ci95"]),
                        "mean_R50_growth_bias_dex": interval_size["mean"],
                        "mean_R50_growth_bias_ci95": json.dumps(interval_size["mean_ci95"]),
                        "median_R50_growth_bias_dex": interval_size["median"],
                        "median_R50_growth_bias_ci95": json.dumps(interval_size["median_ci95"]),
                    })

    # Recover the fair one-output method comparison on the exact 36-tree
    # growth cohort, including the frozen Darwin-calibrated 3D sphere.
    transfer_metrics = []
    transfer_tree_values = {}
    transfer_fields = (
        "light_weighted_IoU",
        "light_completeness",
        "light_precision",
        "log10_mass_ratio_vs_central_bound_dex",
        "log10_R50_ratio_vs_isophote_dex",
        "log10_R50_ratio_vs_central_bound_dex",
    )
    for band in ("V", "g"):
        for method in METHODS:
            for field in transfer_fields:
                tree_values = []
                for rank in RANKS:
                    vals = [float(select_row(rank, 84, band, method, los,
                                             "actual_observer_redshift")[field])
                            for los in ("x", "y", "z")]
                    tree_values.append(float(np.median(vals)))
                stats = bootstrap(np.asarray(tree_values), rng)
                transfer_tree_values[(band, method, field)] = {
                    rank: value for rank, value in zip(RANKS, tree_values)
                }
                transfer_metrics.append({
                    "physical_snapshot": 84,
                    "physical_redshift": SNAP_TO_REDSHIFT[84],
                    "band": band,
                    "method": method,
                    "metric": field,
                    "n_descendant_trees": len(tree_values),
                    "mean": stats["mean"],
                    "mean_ci95": json.dumps(stats["mean_ci95"]),
                    "median": stats["median"],
                    "median_ci95": json.dumps(stats["median_ci95"]),
                })

    transfer_pairs = []
    for band in ("V", "g"):
        for competitor in ("b020_linked_group", "projected_aperture_30_pkpc",
                           "Darwin_calibrated_sphere"):
            for metric in ("light_weighted_IoU", "log10_mass_ratio_vs_central_bound_dex",
                           "log10_R50_ratio_vs_isophote_dex"):
                diffs = np.array([
                    transfer_tree_values[(band, "b010_linked_group", metric)][rank]
                    - transfer_tree_values[(band, competitor, metric)][rank]
                    for rank in RANKS
                ])
                stats = bootstrap(diffs, rng)
                transfer_pairs.append({
                    "band": band,
                    "comparison": f"b010 minus {competitor}",
                    "metric": metric,
                    "n_descendant_trees": len(diffs),
                    "mean_paired_difference": stats["mean"],
                    "mean_paired_difference_ci95": json.dumps(stats["mean_ci95"]),
                    "median_paired_difference": stats["median"],
                    "median_paired_difference_ci95": json.dumps(stats["median_ci95"]),
                })

    paired = []
    for band in ("V", "g"):
        for condition in CONDITIONS:
            b010 = {(r["sample_rank"]): r for r in per_tree
                    if r["band"] == band and r["condition"] == condition
                    and r["method"] == "b010_linked_group"}
            for competitor in ("b020_linked_group", "projected_aperture_30_pkpc",
                               "Darwin_calibrated_sphere", "isophotal_particle_selection"):
                comp = {r["sample_rank"]: r for r in per_tree
                        if r["band"] == band and r["condition"] == condition
                        and r["method"] == competitor}
                mass_diff = np.array([abs(b010[r]["mass_growth_bias_dex"])
                                      - abs(comp[r]["mass_growth_bias_dex"]) for r in RANKS])
                size_diff = np.array([abs(b010[r]["R50_growth_bias_dex"])
                                      - abs(comp[r]["R50_growth_bias_dex"]) for r in RANKS])
                m = bootstrap(mass_diff, rng)
                s = bootstrap(size_diff, rng)
                comp_mass = summaries_by_key[(band, condition, competitor)][0]["mean_absolute_error"]
                b_mass = summaries_by_key[(band, condition, "b010_linked_group")][0]["mean_absolute_error"]
                paired.append({
                    "band": band, "condition": condition,
                    "comparison": f"b010 minus {competitor} in per-tree absolute growth error",
                    "mean_paired_mass_error_difference_dex": m["mean"],
                    "mean_paired_mass_error_difference_ci95": json.dumps(m["mean_ci95"]),
                    "mean_paired_R50_error_difference_dex": s["mean"],
                    "mean_paired_R50_error_difference_ci95": json.dumps(s["mean_ci95"]),
                    "b010_mass_error_reduction_vs_comparator_percent": (
                        100.0 * (1.0 - b_mass / comp_mass) if comp_mass > 0 else float("nan")),
                })

    crossrender = []
    for snap in SNAPS:
        for band in ("V", "g"):
            mass_fraction_diffs = []
            radius_log_diffs = []
            low_detection_by_tree = []
            high_detection_by_tree = []
            total_paired_radius_views = 0
            for rank in RANKS:
                view_mass_diffs = []
                view_radius_diffs = []
                low_detection = []
                high_detection = []
                for los in ("x", "y", "z"):
                    low = select_row(rank, snap, band, "isophotal_particle_selection", los,
                                     "fixed_observer_z0p197")
                    high = [r for r in indexed[(rank, snap, band,
                                                "isophotal_particle_selection", los)]
                            if abs(float(r["observer_redshift"]) - OBSERVER_REDSHIFTS[2]) < 1e-8]
                    if len(high) != 1:
                        raise RuntimeError("Missing high-redshift cross-render row")
                    high = high[0]
                    low_detected = float(low["observed_isophote_area_arcsec2"]) > 0.0
                    high_detected = float(high["observed_isophote_area_arcsec2"]) > 0.0
                    low_detection.append(float(low_detected))
                    high_detection.append(float(high_detected))
                    bound_low = float(low["central_bound_stellar_mass_msun"])
                    bound_high = float(high["central_bound_stellar_mass_msun"])
                    if min(bound_low, bound_high) <= 0:
                        raise RuntimeError("Non-positive central bound mass in cross-render table")
                    low_fraction = float(low["stellar_mass_msun"]) / bound_low
                    high_fraction = float(high["stellar_mass_msun"]) / bound_high
                    view_mass_diffs.append(high_fraction - low_fraction)
                    if low_detected and high_detected:
                        low_r50 = float(low["projected_R50_pkpc"])
                        high_r50 = float(high["projected_R50_pkpc"])
                        if min(low_r50, high_r50) > 0:
                            view_radius_diffs.append(float(np.log10(high_r50 / low_r50)))
                mass_fraction_diffs.append(float(np.median(view_mass_diffs)))
                low_detection_by_tree.append(float(np.mean(low_detection)))
                high_detection_by_tree.append(float(np.mean(high_detection)))
                if view_radius_diffs:
                    radius_log_diffs.append(float(np.median(view_radius_diffs)))
                    total_paired_radius_views += len(view_radius_diffs)

            mass_stat = bootstrap(np.asarray(mass_fraction_diffs), rng)
            radius_stat = (bootstrap(np.asarray(radius_log_diffs), rng, require_full_cohort=False)
                           if radius_log_diffs else None)
            low_det_stat = bootstrap(np.asarray(low_detection_by_tree), rng)
            high_det_stat = bootstrap(np.asarray(high_detection_by_tree), rng)
            crossrender.append({
                "physical_snapshot": snap,
                "physical_redshift": SNAP_TO_REDSHIFT[snap],
                "band": band,
                "contrast": "same particle state at observer z=0.9973 minus observer z=0.1973",
                "low_observer_detection_fraction": low_det_stat["mean"],
                "high_observer_detection_fraction": high_det_stat["mean"],
                "high_observer_detection_fraction_ci95": json.dumps(high_det_stat["mean_ci95"]),
                "mean_change_in_isophotal_mass_fraction_of_bound_central_mass": mass_stat["mean"],
                "mean_isophotal_mass_fraction_change_ci95": json.dumps(mass_stat["mean_ci95"]),
                "n_descendant_trees_mass": len(mass_fraction_diffs),
                "mean_log10_R50_high_over_low_dex_detected_only": (
                    radius_stat["mean"] if radius_stat else ""),
                "R50_log_ratio_ci95_detected_only": (
                    json.dumps(radius_stat["mean_ci95"]) if radius_stat else ""),
                "n_descendant_trees_with_paired_R50": len(radius_log_diffs),
                "n_paired_projection_R50_measurements": total_paired_radius_views,
            })

    write_csv(HERE / "tree_growth_bias_per_tree.csv", per_tree)
    write_csv(HERE / "tree_growth_bias_summary.csv", summary)
    write_csv(HERE / "tree_growth_bias_intervals_summary.csv", interval_summary)
    write_csv(HERE / "tree_growth_bias_paired_comparisons.csv", paired)
    write_csv(HERE / "growth_cohort_z0p2_method_metrics.csv", transfer_metrics)
    write_csv(HERE / "growth_cohort_z0p2_paired_methods.csv", transfer_pairs)
    write_csv(HERE / "crossrender_observational_effects.csv", crossrender)

    # Show mass and size growth errors for actual observer redshift and a fixed
    # low-redshift imaging condition. Every plotted point has one value per tree.
    plt.rcParams.update({"font.size": 10, "axes.labelsize": 10, "xtick.labelsize": 9,
                         "ytick.labelsize": 9, "legend.fontsize": 8,
                         "pdf.fonttype": 42, "ps.fonttype": 42})
    fig, axes = plt.subplots(2, 1, figsize=(7.1, 5.7), sharex=True,
                             gridspec_kw={"hspace": 0.08})
    method_labels = ["Isophote", "$b=0.10$ group", "$b=0.20$ group",
                     "30 pkpc circle", "Darwin sphere"]
    colors = {"V": "#276FBF", "g": "#D95F02"}
    offsets = {("V", CONDITIONS[0]): -0.18, ("g", CONDITIONS[0]): -0.06,
               ("V", CONDITIONS[1]): 0.06, ("g", CONDITIONS[1]): 0.18}
    for ax, stat_kind, ylabel, panel in (
        (axes[0], "mass", r"Mass-growth bias, $\Delta\log_{10}M_{\star}$ (dex)", "(a)"),
        (axes[1], "size", r"Half-light-radius growth bias, $\Delta\log_{10}R_{50}$ (dex)", "(b)"),
    ):
        ax.axhline(0.0, color="0.25", lw=0.9, ls="--", zorder=0)
        ax.text(0.015, 0.93, panel, transform=ax.transAxes, ha="left", va="top",
                fontsize=10, fontweight="bold")
        for i, method in enumerate(METHODS):
            for band in ("V", "g"):
                color = colors[band]
                for condition in CONDITIONS:
                    stat = summaries_by_key[(band, condition, method)][0 if stat_kind == "mass" else 1]
                    y = stat["mean"]
                    lo, hi = stat["mean_ci95"]
                    ax.errorbar(i + offsets[(band, condition)], y,
                                yerr=[[y - lo], [hi - y]], fmt="o", ms=4.5,
                                color=color, mec=color,
                                mfc=(color if condition == CONDITIONS[0] else "white"),
                                capsize=2.2, elinewidth=1.0, capthick=0.8, zorder=3)
        ax.set_ylabel(ylabel)
        ax.grid(axis="y", color="0.88", lw=0.6)
        ax.set_axisbelow(True)
    axes[1].set_xticks(np.arange(len(METHODS)))
    axes[1].set_xticklabels(method_labels)
    axes[1].set_xlabel("Stellar selection method")
    legend = [
        Line2D([0], [0], marker="o", color="none", markerfacecolor=colors["V"],
               markeredgecolor=colors["V"], label="$V$, actual observer redshift", markersize=5),
        Line2D([0], [0], marker="o", color="none", markerfacecolor=colors["g"],
               markeredgecolor=colors["g"], label="$g$, actual observer redshift", markersize=5),
        Line2D([0], [0], marker="o", color="0.2", markerfacecolor="white",
               markeredgecolor="0.2", label="Fixed observer redshift $z=0.1973$", markersize=5),
    ]
    axes[0].legend(handles=legend, loc="best", frameon=False, ncol=1)
    fig.savefig(FIGURE, bbox_inches="tight")
    plt.close(fig)

    result = {
        "simulation": "TNG300-1",
        "n_independent_descendant_trees": len(RANKS),
        "statistical_unit": "one descendant tree; median across three projections before bootstrap",
        "snapshots": {str(s): {"redshift": SNAP_TO_REDSHIFT[s]} for s in SNAPS},
        "observer_conditions": {
            "actual_observer_redshift": "observer redshift equals the physical snapshot redshift",
            "fixed_observer_z0p197": "all physical states rendered at observer redshift 0.1973",
        },
        "crossrender_effect": "same physical particle distribution rendered at observer z=0.1973 and 0.9973",
        "mass_definition": "sum of current stellar particle masses in selected particle set",
        "reference_definition": "SubLink central subhalo star ParticleIDs matched exactly to parent FoF cutout",
        "maximum_absolute_log10_catalog_to_particle_mass_ratio_dex": float(max_catalog_mass_log_error),
        "bootstrap_seed": SEED,
        "bootstrap_replicates": N_BOOT,
        "summary_csv": "tree_growth_bias_summary.csv",
        "paired_comparisons_csv": "tree_growth_bias_paired_comparisons.csv",
        "crossrender_csv": "crossrender_observational_effects.csv",
        "figure": str(FIGURE.relative_to(PAPER)),
        "subLink_reference": "Rodriguez-Gomez et al. 2015, MNRAS, 449, 49",
    }
    summary_json = HERE / "tree_growth_bias_analysis.json"
    if summary_json.exists():
        raise FileExistsError(f"Refusing to overwrite {summary_json.name}")
    summary_json.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({"n_trees": len(RANKS),
                      "max_mass_audit_error_dex": max_catalog_mass_log_error,
                      "figure": str(FIGURE)}, indent=2), flush=True)


if __name__ == "__main__":
    main()
