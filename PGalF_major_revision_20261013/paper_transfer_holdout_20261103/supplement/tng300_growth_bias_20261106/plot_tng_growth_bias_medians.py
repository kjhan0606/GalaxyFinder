#!/usr/bin/env python3
"""Plot median TNG300-1 growth biases and descendant-tree bootstrap intervals."""
from __future__ import annotations

import csv
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np

HERE = Path(__file__).resolve().parent
PAPER = HERE.parents[1]
INPUT = HERE / "tree_growth_bias_summary.csv"
OUTPUT = PAPER / "figures/tng300_tree_growth_bias_median.pdf"


def main() -> None:
    with INPUT.open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    methods = [
        "isophotal_particle_selection",
        "b010_linked_group",
        "b020_linked_group",
        "projected_aperture_30_pkpc",
        "Darwin_calibrated_sphere",
    ]
    conditions = ("actual_observer_redshift", "fixed_observer_z0p197")
    by_key = {(r["band"], r["condition"], r["method"]): r for r in rows}
    colors = {"V": "#276FBF", "g": "#D95F02"}
    offsets = {
        ("V", conditions[0]): -0.18,
        ("g", conditions[0]): -0.06,
        ("V", conditions[1]): 0.06,
        ("g", conditions[1]): 0.18,
    }

    plt.rcParams.update({
        "font.size": 10,
        "axes.labelsize": 10,
        "xtick.labelsize": 9,
        "ytick.labelsize": 9,
        "legend.fontsize": 10,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    })
    fig, axes = plt.subplots(
        2, 1, figsize=(7.1, 5.7), sharex=True,
        gridspec_kw={"hspace": 0.17},
    )
    panel_specs = (
        ("median_mass_growth_bias_dex", "median_mass_growth_bias_ci95",
         "Mass-growth bias (dex)", "(a)"),
        ("median_R50_growth_bias_dex", "median_R50_growth_bias_ci95",
         r"Projected $R_{50}$ growth bias (dex)", "(b)"),
    )
    for ax, (value_key, interval_key, ylabel, panel) in zip(axes, panel_specs):
        ax.axhline(0.0, color="0.25", lw=0.9, ls="--", zorder=0)
        ax.text(0.015, 0.95, panel, transform=ax.transAxes, ha="left", va="top",
                fontsize=10, fontweight="bold")
        for i, method in enumerate(methods):
            for band in ("V", "g"):
                for condition in conditions:
                    record = by_key[(band, condition, method)]
                    value = float(record[value_key])
                    lo, hi = json.loads(record[interval_key])
                    color = colors[band]
                    ax.errorbar(
                        i + offsets[(band, condition)], value,
                        yerr=[[value - lo], [hi - value]], fmt="o", ms=4.5,
                        color=color, mec=color,
                        mfc=color if condition == conditions[0] else "white",
                        capsize=2.2, elinewidth=1.0, capthick=0.8, zorder=3,
                    )
        ax.set_ylabel(ylabel, labelpad=8)
        ax.grid(axis="y", color="0.88", lw=0.6)
        ax.set_axisbelow(True)

    method_labels = [
        "Isophote", "$b=0.10$ group", "$b=0.20$ group",
        "30 pkpc circle", "Darwin sphere",
    ]
    axes[1].set_xticks(np.arange(len(methods)))
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
    axes[0].legend(handles=legend, loc="upper right", frameon=False, ncol=1,
                   borderaxespad=0.35, handletextpad=0.45, labelspacing=0.25)
    fig.subplots_adjust(left=0.16, right=0.99, top=0.99, bottom=0.13)
    fig.savefig(OUTPUT, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
