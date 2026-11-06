#!/usr/bin/env python3
"""Plot branch-level stellar-mass growth errors and mass recovery histories."""
from __future__ import annotations

import csv
import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D

HERE = Path(__file__).resolve().parent
PAPER = HERE.parents[1]
TNG = PAPER / "supplement" / "tng300_growth_bias_20261106"
AUDIT = TNG / "growth_r1r4_audit_20261106"
OUTPUT = HERE / "branch_level_growth_recovery.pdf"
METHODS = {
    "isophotal_particle_selection": "Isophote",
    "b012_linked_group": r"$b_{\mathrm{link}}=0.12$",
    "b010_linked_group": r"$b_{\mathrm{link}}=0.10$",
    "b020_linked_group": r"$b_{\mathrm{link}}=0.20$",
}
COLORS = {
    "isophotal_particle_selection": "#C65D3A",
    "b012_linked_group": "#276FBF",
    "b010_linked_group": "#276FBF",
    "b020_linked_group": "#628A48",
}
ENDPOINTS = ((50, r"$z\simeq1.0$"), (67, r"$z\simeq0.5$"), (84, r"$z\simeq0.2$"))


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as fh:
        return list(csv.DictReader(fh))


def scatter_growth(ax, rows, band, calibrated_method, ylabel=True):
    selected_methods = ("isophotal_particle_selection", calibrated_method, "b020_linked_group")
    rng = np.random.default_rng(20261106 + (0 if band == "V" else 1))
    x = np.arange(3)
    for i, method in enumerate(selected_methods):
        q = [r for r in rows if r["band"] == band
             and r["condition"] == "actual_observer_redshift"
             and r["method"] == method]
        q.sort(key=lambda r: int(r["sample_rank"]))
        vals = np.asarray([float(r["growth_bias_stellar_mass_dex"]) for r in q])
        detected = np.asarray([bool(int(r["detected_at_both_endpoints_in_all_three_views"]))
                               for r in q])
        jitter = rng.uniform(-0.105, 0.105, len(q))
        color = COLORS[method]
        ax.scatter(x[i] + jitter[detected], vals[detected], s=20, color=color,
                   edgecolor="white", linewidth=0.35, alpha=0.85, zorder=3)
        if np.any(~detected):
            ax.scatter(x[i] + jitter[~detected], vals[~detected], s=31, marker="^",
                       facecolor="none", edgecolor=color, linewidth=0.9, zorder=4)
        med = float(np.median(vals))
        ax.plot([x[i] - 0.19, x[i] + 0.19], [med, med], color="0.1", lw=1.2, zorder=5)
    ax.axhline(0.0, color="0.4", lw=0.8, ls="--", zorder=0)
    ax.set_xticks(x, [METHODS[m] for m in selected_methods])
    ax.tick_params(axis="x", labelrotation=12)
    ax.grid(axis="y", color="0.9", linewidth=0.6)
    ax.set_axisbelow(True)
    if ylabel:
        ax.set_ylabel(r"$\Delta\log_{10} M_{\star}$ (dex)")
    title_position = (0.97, 0.80) if band == "V" else (0.35, 0.48)
    title_align = "right" if band == "V" else "left"
    ax.text(*title_position, "V" if band == "V" else "SDSS $g$",
            transform=ax.transAxes, ha=title_align, va="top", fontsize=9.5,
            fontweight="bold", bbox=dict(facecolor="white", edgecolor="none", alpha=.9, pad=1.2))


def add_g_zoom(ax, rows):
    inset = ax.inset_axes([0.39, 0.52, 0.57, 0.43])
    methods = ("isophotal_particle_selection", "b010_linked_group", "b020_linked_group")
    rng = np.random.default_rng(20261109)
    for i, method in enumerate(methods):
        q = [r for r in rows if r["band"] == "g"
             and r["condition"] == "actual_observer_redshift" and r["method"] == method]
        q.sort(key=lambda r: int(r["sample_rank"]))
        vals = np.asarray([float(r["growth_bias_stellar_mass_dex"]) for r in q])
        detected = np.asarray([bool(int(r["detected_at_both_endpoints_in_all_three_views"]))
                               for r in q])
        jitter = rng.uniform(-0.105, 0.105, len(q))
        inset.scatter(i + jitter[detected], vals[detected], s=10, color=COLORS[method],
                      edgecolor="white", linewidth=0.25, alpha=0.8)
        if np.any(~detected):
            inset.scatter(i + jitter[~detected], vals[~detected], s=15, marker="^",
                          facecolor="none", edgecolor=COLORS[method], linewidth=0.7)
        inset.plot([i - .18, i + .18], [np.median(vals), np.median(vals)], color="0.1", lw=0.8)
    inset.axhline(0, color="0.4", lw=0.6, ls="--")
    inset.set_xlim(-.45, 2.45)
    inset.set_ylim(-.16, .22)
    inset.set_xticks([0, 1, 2], ["Isophote", "$b=0.10$", "$b=0.20$"])
    inset.tick_params(axis="x", labelsize=6.2, rotation=10, pad=0,
                      top=True, bottom=True, direction="in")
    inset.tick_params(axis="y", labelsize=6.4, left=True, right=True,
                      direction="in")
    inset.text(.97, .92, "(b)", transform=inset.transAxes, ha="right", va="top",
               fontsize=7.2, fontweight="bold",
               bbox=dict(facecolor="white", edgecolor="none", alpha=.9, pad=.4))
    inset.grid(axis="y", color="0.9", linewidth=0.45)
    inset.set_facecolor("white")


def recovery_panel(ax, rows, band, method_names):
    for method in method_names:
        per_snap = []
        q16, medians, q84 = [], [], []
        for snap, _label in ENDPOINTS:
            vals = np.asarray([
                float(r["median_selected_to_central_mass_ratio"])
                for r in rows if r["band"] == band
                and r["observer_condition"] == "actual_observer_redshift"
                and r["snapshot"] == str(snap) and r["method"] == method
            ])
            q16.append(float(np.quantile(vals, .16)))
            medians.append(float(np.median(vals)))
            q84.append(float(np.quantile(vals, .84)))
        x = np.arange(3)
        color = COLORS[method]
        ax.plot(x, medians, marker="o", ms=3.6, lw=1.25, color=color,
                label=METHODS[method])
        ax.fill_between(x, q16, q84, color=color, alpha=.13, linewidth=0)
    ax.set_xticks(np.arange(3), [label for _, label in ENDPOINTS])
    ax.set_ylim(0.0, 1.0)
    ax.set_ylabel(r"$M_{\rm selected}/M_{\rm central}$")
    ax.grid(axis="y", color="0.9", linewidth=0.6)
    ax.set_axisbelow(True)
    ax.text(.025, .97, "V" if band == "V" else "SDSS $g$",
            transform=ax.transAxes, ha="left", va="top", fontsize=9.5,
            fontweight="bold", bbox=dict(facecolor="white", edgecolor="none", alpha=.9, pad=1.2))
    ax.legend(loc="lower right", frameon=True, framealpha=.92, fontsize=6.8,
              borderpad=.25, labelspacing=.2, handlelength=1.2)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--force", action="store_true",
                        help="replace the existing generated figure")
    args = parser.parse_args()
    if OUTPUT.exists() and not args.force:
        raise FileExistsError(f"Refusing to overwrite {OUTPUT}")
    growth = read_csv(AUDIT / "r1r3_growth_bias_per_tree.csv")
    recovery = read_csv(AUDIT / "r2_mass_recovery_per_tree.csv")
    plt.rcParams.update({"font.size": 8, "axes.labelsize": 8.5,
                         "xtick.labelsize": 7.2, "ytick.labelsize": 7.5,
                         "xtick.direction": "in", "ytick.direction": "in",
                         "pdf.fonttype": 42, "ps.fonttype": 42})
    fig, axes = plt.subplots(2, 2, figsize=(7.05, 5.3),
                             gridspec_kw={"hspace": .28, "wspace": .23})
    scatter_growth(axes[0, 0], growth, "V", "b012_linked_group")
    scatter_growth(axes[0, 1], growth, "g", "b010_linked_group", ylabel=False)
    axes[0, 0].set_ylim(-.18, .43)
    axes[0, 1].set_ylim(-.25, 3.15)
    add_g_zoom(axes[0, 1], growth)
    handles = [
        Line2D([], [], marker="o", ls="none", color="0.25", markerfacecolor="0.25",
               markeredgecolor="white", label="Central pixel detected"),
        Line2D([], [], marker="^", ls="none", color="0.25", markerfacecolor="none",
               markeredgecolor="0.25", label="Fallback in an endpoint view"),
    ]
    recovery_panel(axes[1, 0], recovery, "V",
                   ("isophotal_particle_selection", "b012_linked_group", "b020_linked_group"))
    recovery_panel(axes[1, 1], recovery, "g",
                   ("isophotal_particle_selection", "b010_linked_group", "b020_linked_group"))
    for ax in axes.flat:
        ax.tick_params(top=True, bottom=True, left=True, right=True, direction="in")
    axes[1, 0].set_xlabel("Simulation output")
    axes[1, 1].set_xlabel("Simulation output")
    for ax, label in zip(axes.flat, ("(a)", "(b)", "(c)", "(d)")):
        if label == "(b)":
            continue
        x = .975
        ax.text(x, .97, label, transform=ax.transAxes, ha="right", va="top",
                fontsize=8.5, fontweight="bold",
                bbox=dict(facecolor="white", edgecolor="none", alpha=.9, pad=.5))
    fig.subplots_adjust(left=.105, right=.99, top=.985, bottom=.12)
    fig.savefig(OUTPUT, bbox_inches="tight")
    plt.close(fig)
    print(OUTPUT)


if __name__ == "__main__":
    main()
