#!/usr/bin/env python3
"""Plot updated TNG growth biases including the calibrated V-band b=0.12."""
from __future__ import annotations

import csv
import json
import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np

HERE = Path(__file__).resolve().parent
PAPER = HERE.parents[1]
DEFAULT_SUMMARY = HERE / "growth_r1r4_audit_20261106" / "r1r3_growth_bias_summary.csv"
DEFAULT_OUTPUT = PAPER / "figures/tng300_tree_growth_bias_r1r4d.pdf"
DETECTED_COHORT = "central_pixel_detected_all_views_at_both_endpoints"
METHODS = (
    "isophotal_particle_selection",
    "b010_linked_group",
    "b012_linked_group",
    "b020_linked_group",
    "projected_aperture_30_pkpc",
    "Darwin_calibrated_sphere",
)
CONDITIONS = ("actual_observer_redshift", "fixed_observer_z0p197")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--summary", type=Path, default=DEFAULT_SUMMARY)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--force", action="store_true",
                        help="replace the existing generated figure")
    args = parser.parse_args()
    output = args.output.resolve()
    if output.exists() and not args.force:
        raise FileExistsError(f"Refusing to overwrite {output}")
    with args.summary.open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    indexed = {(row["cohort"], row["band"], row["observer_condition"],
                row["method"], row["measurement"]): row for row in rows}
    colors = {"V": "#276FBF", "g": "#D95F02"}
    offsets = {
        ("V", CONDITIONS[0]): -0.24,
        ("g", CONDITIONS[0]): -0.08,
        ("V", CONDITIONS[1]): 0.08,
        ("g", CONDITIONS[1]): 0.24,
    }
    plt.rcParams.update({
        "font.size": 9,
        "axes.labelsize": 9,
        "xtick.labelsize": 8,
        "ytick.labelsize": 9,
        "legend.fontsize": 9,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    })
    fig, axes = plt.subplots(2, 1, figsize=(7.1, 5.9), sharex=True,
                             gridspec_kw={"hspace": 0.0})
    panel_specs = (
        ("stellar_mass", r"$\Delta\log_{10}M_{\star}$ (dex)", "(a)"),
        ("half_light_radius", r"$\Delta\log_{10}R_{50}$ (dex)", "(b)"),
    )
    for ax, (measurement, ylabel, panel) in zip(axes, panel_specs):
        ax.axhline(0.0, color="0.25", lw=0.9, ls="--", zorder=0)
        title_x, title_ha = (0.985, "right") if panel == "(b)" else (0.015, "left")
        ax.text(title_x, 0.95, panel, transform=ax.transAxes, ha=title_ha, va="top",
                fontsize=9, fontweight="bold")
        for i, method in enumerate(METHODS):
            for band in ("V", "g"):
                for condition in CONDITIONS:
                    row = indexed[("all_36_branches", band, condition,
                                   method, measurement)]
                    value = float(row["median_bias_dex"])
                    lo, hi = json.loads(row["median_bias_ci95"])
                    color = colors[band]
                    ax.errorbar(
                        i + offsets[(band, condition)], value,
                        yerr=[[value - lo], [hi - value]], fmt="o", ms=4.5,
                        color=color, mec=color,
                        mfc=color if condition == CONDITIONS[0] else "white",
                        capsize=2.2, elinewidth=1.0, capthick=0.8, zorder=3,
                    )
            if method == METHODS[0]:
                row = indexed[(DETECTED_COHORT, "g", CONDITIONS[0],
                               method, measurement)]
                value = float(row["median_bias_dex"])
                lo, hi = json.loads(row["median_bias_ci95"])
                ax.errorbar(
                    i + 0.62, value,
                    yerr=[[value - lo], [hi - value]], fmt="D", ms=4.4,
                    color=colors["g"], mec="0.15", mfc=colors["g"],
                    capsize=2.2, elinewidth=1.0, capthick=0.8, zorder=4,
                )
        ax.set_ylabel(ylabel, labelpad=8)
        ax.tick_params(axis="y", left=True, right=True, direction="in")
        ax.grid(axis="y", color="0.88", lw=0.6)
        ax.set_axisbelow(True)

    labels = ("Isophote", r"$b_{\mathrm{link}}=0.10$", r"$b_{\mathrm{link}}=0.12$", r"$b_{\mathrm{link}}=0.20$",
              "30 pkpc", "80/60 pkpc sphere")
    axes[1].set_xticks(np.arange(len(METHODS)))
    axes[1].set_xticklabels(labels, rotation=18, ha="right", rotation_mode="anchor")
    axes[1].set_xlabel("Method")
    axes[0].tick_params(axis="x", bottom=True, top=True, labelbottom=False, labeltop=False, direction="in")
    axes[1].tick_params(axis="x", bottom=True, top=True, labeltop=False, direction="in")
    legend = [
        Line2D([0], [0], marker="o", color="none", markerfacecolor=colors["V"],
               markeredgecolor=colors["V"], label="$V$, physical $z_{\\rm obs}$", markersize=5),
        Line2D([0], [0], marker="o", color="none", markerfacecolor=colors["g"],
               markeredgecolor=colors["g"], label="$g$, physical $z_{\\rm obs}$", markersize=5),
        Line2D([0], [0], marker="D", color="none", markerfacecolor=colors["g"],
               markeredgecolor="0.15",
               label="$g$, central detections (29 branches)", markersize=5),
        Line2D([0], [0], marker="o", color="0.2", markerfacecolor="white",
               markeredgecolor="0.2", label="Fixed $z_{\\rm obs}=0.1973$", markersize=5),
    ]
    axes[0].legend(handles=legend, loc="upper right", frameon=False, ncol=1,
                   borderaxespad=0.35, handletextpad=0.45, labelspacing=0.25)
    fig.subplots_adjust(left=0.16, right=0.99, top=0.99, bottom=0.17)
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, bbox_inches="tight")
    plt.close(fig)
    print(output)


if __name__ == "__main__":
    main()
