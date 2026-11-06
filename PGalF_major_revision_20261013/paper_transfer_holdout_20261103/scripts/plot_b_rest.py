#!/usr/bin/env python3
"""Regenerate the archived rest-frame overlap figure with inward ticks."""
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "figures" / "b_rest.pdf"
Z = (0.190, 0.227, 0.500, 1.002, 1.990, 3.002)
COLORS = dict(zip(Z, ("#1F4E79", "#148F9A", "#5B8C3E", "#E0A100", "#D4652F", "#8E3A73")))
STYLES = ("-", "--", "-.", ":")


def main():
    rest = defaultdict(list)
    masses = {}
    for line in (ROOT / "supplement" / "rest_frame.tsv").read_text().splitlines():
        if line.startswith("#") or not line.strip():
            continue
        p = line.split()
        if p[1] == "P":
            masses[(p[0], int(p[3]))] = float(p[7])
        if p[1] == "I":
            rest[(p[0], int(p[3]), p[5], round(float(p[2]), 3))].append(
                [p[4], float(p[6]), float(p[7])]
            )

    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 8,
                         "xtick.direction": "in", "ytick.direction": "in",
                         "pdf.fonttype": 42})
    fig, axes = plt.subplots(2, 1, figsize=(3.45, 4.5), sharex=True, sharey=True)
    for index, band in enumerate(("V", "g")):
        ax = axes[index]
        ax.tick_params(top=True, bottom=True, left=True, right=True, direction="in")
        band_items = sorted((k, v) for k, v in rest.items() if k[2] == band)
        for series_index, ((step, nstar, _band, redshift), rows) in enumerate(band_items):
            rows = sorted(rows, key=lambda row: row[1])
            x = np.asarray(sorted({row[1] for row in rows}))
            curves = np.asarray([[row[2] for row in rows if row[0] == los]
                                 for los in "xyz"])
            color = COLORS[redshift]
            mass = masses[(step, nstar)]
            ax.plot(x, np.median(curves, axis=0), color=color,
                    ls=STYLES[series_index % len(STYLES)], lw=1.4,
                    label=fr"$z={redshift:.3f}$, $M={mass/1e13:.3f}" +
                          r"\times10^{13}$")
            ax.fill_between(x, curves.min(axis=0), curves.max(axis=0),
                            color=color, alpha=.09)
        ax.set(xlim=(.01, .51), ylim=(0, 1.035), ylabel=r"$\mathrm{IoU}$")
        ax.grid(alpha=.18)
        ax.text(.03, .96, f"({chr(97 + index)})", transform=ax.transAxes,
                va="top", fontsize=9)
        ax.text(.97, .65, "$V$" if band == "V" else "SDSS $g$",
                transform=ax.transAxes, ha="right", va="top", fontsize=9)
    axes[1].set_xlabel(r"$b_{\mathrm{link}}$")
    axes[1].legend(loc="lower right", fontsize=9, frameon=False)
    fig.subplots_adjust(left=.17, right=.99, bottom=.10, top=.99, hspace=0)
    # Keep the inward tick at the join; show the shared-edge number once.
    axes[0].yaxis.get_major_ticks()[0].label1.set_visible(False)
    axes[0].tick_params(labelbottom=False)
    fig.savefig(OUT, bbox_inches="tight")
    plt.close(fig)
    print(OUT)


if __name__ == "__main__":
    main()
