#!/usr/bin/env python3
"""Recreate the fixed-structure CAMELS redshift-transfer figure."""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
DATA = HERE / "condition_metrics.csv"
OUTPUT = HERE / "camels_multiz_transfer_summary_v2.pdf"
REDSHIFTS = (0.190438, 0.500319, 1.002070, 1.990001)
METHODS = (
    ("FoF", 0.10, "#2474B5", "o", r"$b_{\rm link}=0.10$"),
    ("FoF", 0.20, "#D95F02", "s", r"$b_{\rm link}=0.20$"),
    ("aperture30pkpc", None, "#159C7A", "^", "30 pkpc aperture"),
)


def at_redshift(df, z):
    return df[np.isclose(df["observer_redshift"], z, rtol=0, atol=2e-5)]


def main():
    df = pd.read_csv(DATA)
    df = df[df["threshold_case"] == "base"]
    plt.rcParams.update({"font.size": 8, "axes.labelsize": 9,
                         "xtick.direction": "in", "ytick.direction": "in",
                         "pdf.fonttype": 42})
    fig, axes = plt.subplots(2, 2, figsize=(7.05, 4.7), sharex="col", sharey="row")
    x = np.arange(len(REDSHIFTS))
    for col, band in enumerate(("V", "g")):
        ax = axes[0, col]
        for method, b, color, marker, label in METHODS:
            rows = []
            for z in REDSHIFTS:
                q = at_redshift(df, z)
                q = q[(q["band"] == band) & (q["method"] == method)]
                if b is not None:
                    q = q[np.isclose(pd.to_numeric(q["b_star"], errors="coerce"), b)]
                if len(q) != 1:
                    raise ValueError(f"Expected one row for {band}, {method}, b={b}, z={z}")
                rows.append(q.iloc[0])
            y = np.asarray([r["mean_host_median_iou"] for r in rows], dtype=float)
            lo = np.asarray([r["mean_host_median_iou_box_then_host_bootstrap_95_low"]
                             for r in rows], dtype=float)
            hi = np.asarray([r["mean_host_median_iou_box_then_host_bootstrap_95_high"]
                             for r in rows], dtype=float)
            ax.errorbar(x, y, yerr=np.vstack((y - lo, hi - y)), marker=marker,
                        ms=4, lw=1.25, capsize=2, color=color, label=label)

        fail_rows = []
        for z in REDSHIFTS:
            q = at_redshift(df, z)
            q = q[q["band"] == band]
            q = q[(q["method"] == "FoF") &
                  np.isclose(pd.to_numeric(q["b_star"], errors="coerce"), .10)]
            fail_rows.append(q.iloc[0])
        y = np.asarray([r["fraction_halos_with_any_failed_projection"] for r in fail_rows])
        axes[1, col].plot(x, y, marker="D", ms=3.5, lw=1.1, ls="--",
                          color="0.3", label="At least one failed view")
        for row in range(2):
            a = axes[row, col]
            a.set_xticks(x, [f"{z:.1f}" for z in REDSHIFTS])
            a.set_ylim(0, 1)
            a.grid(axis="y", alpha=.18)
            a.tick_params(top=True, bottom=True, left=True, right=True, direction="in",
                          labelbottom=(row == 1), labeltop=False,
                          labelleft=(col == 0), labelright=False)
            a.text(.025, .96, f"({chr(97 + row*2 + col)}) {band} band",
                   transform=a.transAxes, ha="left", va="top", fontsize=9)
        if col == 0:
            axes[0, col].set_ylabel("Mean host overlap")
            axes[1, col].set_ylabel("Hosts with a failed view")
        axes[1, col].set_xlabel(r"Observer redshift $z_{\rm obs}$")

    # Preserve both panels' ticks at the shared edge, but print the boundary
    # value only once instead of overlapping 0.0 and 1.0 at the join.
    axes[0, 0].yaxis.get_major_ticks()[0].label1.set_visible(False)
    axes[0, 0].legend(loc="upper right", frameon=True, framealpha=.9, fontsize=9)
    fig.subplots_adjust(left=.09, right=.995, bottom=.12, top=.99,
                        wspace=0, hspace=0)
    fig.savefig(OUTPUT, bbox_inches="tight")
    plt.close(fig)
    print(OUTPUT)


if __name__ == "__main__":
    main()
