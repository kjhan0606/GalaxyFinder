#!/usr/bin/env python3
"""Regenerate the cluster illustration from saved display arrays.

No particle catalogue or linking library is needed. Unused full-resolution
stellar arrays are omitted from the compact display archive.
"""
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse, Rectangle
from matplotlib.ticker import MaxNLocator
from paper_figure_style import save_paper_figure
ROOT = Path(__file__).resolve().parents[1]
B_SHOW = (0.05, 0.08, 0.10, 0.15, 0.20)

def render(path):
    """Draw from the saved arrays. Safe to repeat without the halo read."""
    d = np.load(path)
    half = float(d["half"])
    wide = float(d["wide"])
    as_per_pkpc = float(d["as_per_pkpc"])
    xg = d["xg"]
    region = d["region"].astype(float)
    plt.rcParams.update({
        "font.family": "DejaVu Sans",
        "font.size": 8,
        "axes.labelsize": 8,
        "axes.titlesize": 8,
        "xtick.labelsize": 7,
        "ytick.labelsize": 7,
    })
    # Both rows share edges. Panel (f) sits against (e) and under (c).
    fig_w, fig_h = 6.83, 4.68
    pw, ph = 2.02, 2.02
    left, bottom = 0.62, 0.42
    hgap, vgap = 0.0, 0.0
    fig = plt.figure(figsize=(fig_w, fig_h))
    xs = [
        left, left + pw, left + 2.0 * pw,
        left, left + pw, left + 2.0 * pw + hgap,
    ]
    ys = [bottom + ph + vgap] * 3 + [bottom] * 3

    def place(i):
        return fig.add_axes([
            xs[i] / fig_w, ys[i] / fig_h, pw / fig_w, ph / fig_h,
        ])

    axes = [place(i) for i in range(6)]
    # prune drops the tick label that would collide across a shared edge
    x_prune = ["upper", "both", "lower", "upper", "lower", None]
    for ip, bpar in enumerate(B_SHOW):
        ax = axes[ip]
        ax.scatter(d["xout%d" % ip], d["yout%d" % ip], s=1.2, c="#c47b4a",
                   alpha=0.28, linewidths=0, rasterized=True,
                   label="outside the group")
        ax.scatter(d["xin%d" % ip], d["yin%d" % ip], s=1.2, c="#1f4e79",
                   alpha=0.45, linewidths=0, rasterized=True,
                   label="central group")
        ax.contour(xg, xg, region, levels=[0.5], colors="#111111", linewidths=0.7)
        ax.set_xlim(-half, half)
        ax.set_ylim(-half, half)
        ax.set_aspect("equal")
        ax.tick_params(top=False, right=False, direction="in")
        ax.xaxis.set_major_locator(MaxNLocator(nbins=4, prune=x_prune[ip], steps=[1, 2, 5, 10]))
        ax.yaxis.set_major_locator(MaxNLocator(nbins=4, prune="both", steps=[1, 2, 5, 10]))
        if ip % 3 != 0:
            ax.tick_params(labelleft=False)
        else:
            ax.set_ylabel("arcsec")
        if ip < 2:
            ax.spines["bottom"].set_visible(False)
            ax.tick_params(bottom=False, labelbottom=False)
        elif ip == 2:
            ax.tick_params(bottom=False, labelbottom=False)
        elif ip < 5:
            ax.set_xlabel("arcsec")
        ax.text(
            0.04, 0.96, rf"$b={bpar:.2f}$", transform=ax.transAxes,
            ha="left", va="top", fontsize=8, color="black",
            bbox=dict(boxstyle="square,pad=0.12", facecolor="white",
                      edgecolor="none", alpha=0.88),
        )
        ax.text(
            0.96, 0.96, f"({chr(ord('a') + ip)})", transform=ax.transAxes,
            ha="right", va="top", fontsize=8, color="black",
            bbox=dict(boxstyle="square,pad=0.08", facecolor="white",
                      edgecolor="none", alpha=0.88),
        )
        if ip == 0:
            ax.legend(loc="lower left", frameon=False, fontsize=6.5,
                      markerscale=4, handletextpad=0.2)
    # 100 proper kpc. At this redshift the bar is about 30 arcsec.
    bar_as = 100.0 * as_per_pkpc
    ax0 = axes[0]
    x1 = half - 48.0
    x0 = x1 - bar_as
    y0 = -half + 22.0
    ax0.plot([x0, x1], [y0, y0], color="black", lw=1.1, solid_capstyle="butt")
    ax0.text((x0 + x1) / 2.0, y0 + 5.0, "100 kpc", ha="center", va="bottom",
             fontsize=6.5)

    ax = axes[5]
    hist = np.ma.masked_less_equal(d["hist"].T, 0.0)
    show = np.ma.log10(hist)
    # log10 of stellar mass per bin. The floor near 4.73 is one particle
    # and would paint the field one colour; clipping near 9 flattens
    # every core to the same yellow. These limits keep filaments, galaxy
    # bodies, and cores on separate parts of the scale.
    vmin = 5.55
    vmax = 10.45
    cmap = plt.get_cmap("afmhot").copy()
    cmap.set_bad("#000000")
    ax.set_facecolor("#000000")
    mass_image = ax.imshow(show, origin="lower", extent=[-wide, wide, -wide, wide],
              cmap=cmap, vmin=vmin, vmax=vmax,
              interpolation="nearest", aspect="equal", zorder=1)
    cax = fig.add_axes([.70, .050, .25, .016])
    cb = fig.colorbar(mass_image, cax=cax, orientation="horizontal")
    cb.set_ticks([6, 8, 10]); cb.ax.tick_params(labelsize=6, length=2, pad=1)
    cb.set_label(r"$\log_{10}(M_\star/[h^{-1}M_\odot])$ per bin", fontsize=6, labelpad=1)
    scale = float(d["scale"])
    mean = d["ell_mean"]
    cx = mean[:, 0] * scale
    cy = mean[:, 1] * scale
    bcg = d["ell_bcg"]
    # Median semi-major axis is 0.89 arcsec, below a hairline on this
    # 1 Mpc field. Draw five times the second-moment semi-axes so the
    # size ranking stays visible. Axis ratio and angle are unchanged.
    show_ell = 5.0
    order = np.argsort(np.asarray(d["ell_a"], float))
    for i in order:
        if bool(bcg[i]):
            continue
        semi_a = float(d["ell_a"][i]) * scale * show_ell
        semi_b = float(d["ell_b"][i]) * scale * show_ell
        if not np.isfinite(semi_a) or semi_a <= 0.0:
            continue
        if not np.isfinite(semi_b) or semi_b <= 0.0:
            continue
        ax.add_patch(Ellipse(
            (cx[i], cy[i]), 2.0 * semi_a, 2.0 * semi_b,
            angle=float(d["ell_ang"][i]), facecolor="none",
            edgecolor="#d7f6ff", linewidth=0.28, zorder=4,
        ))
    i_bcg = int(np.flatnonzero(bcg)[0])
    bcg_a = float(d["ell_a"][i_bcg]) * scale * show_ell
    bcg_b = float(d["ell_b"][i_bcg]) * scale * show_ell
    bcg_ang = float(d["ell_ang"][i_bcg])
    for edge, width, zo in (("#000000", 1.35, 8), ("#ffffff", 0.7, 9)):
        ax.add_patch(Ellipse(
            (cx[i_bcg], cy[i_bcg]), 2.0 * bcg_a, 2.0 * bcg_b, angle=bcg_ang,
            facecolor="none", edgecolor=edge, linewidth=width, zorder=zo,
        ))
    ax.add_patch(Rectangle(
        (-half, -half), 2 * half, 2 * half, fill=False,
        edgecolor="#000000", linewidth=1.15, linestyle=(0, (2.4, 1.5)), zorder=6,
    ))
    ax.add_patch(Rectangle(
        (-half, -half), 2 * half, 2 * half, fill=False,
        edgecolor="#ffffff", linewidth=0.6, linestyle=(0, (2.4, 1.5)), zorder=7,
    ))
    ax.set_xlim(-wide, wide)
    ax.set_ylim(-wide, wide)
    ax.set_aspect("equal")
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_color("#000000")
        spine.set_linewidth(0.6)
    ax.text(
        0.96, 0.96, "(f)", transform=ax.transAxes,
        ha="right", va="top", fontsize=8, color="white",
        bbox=dict(boxstyle="square,pad=0.08", facecolor="black",
                  edgecolor="none", alpha=0.55),
        zorder=8,
    )
    # 1 proper Mpc. A 100 kpc bar is only ~30 arcsec on this field.
    wide_bar = 1000.0 * as_per_pkpc
    x0 = -wide + 0.10 * wide
    yb = -wide + 0.12 * wide
    x1 = x0 + wide_bar
    ax.plot([x0, x1], [yb, yb], color="black", lw=2.0, solid_capstyle="butt",
            zorder=8)
    ax.plot([x0, x1], [yb, yb], color="white", lw=1.15, solid_capstyle="butt",
            zorder=9)
    ax.text(
        (x0 + x1) / 2.0, yb + 0.035 * wide, "1 Mpc",
        ha="center", va="bottom", fontsize=6.5, color="white", zorder=9,
        bbox=dict(boxstyle="square,pad=0.08", facecolor="black",
                  edgecolor="none", alpha=0.55),
    )
    print("  bar_arcsec", bar_as, "wide_bar_arcsec", wide_bar,
          "n_ell", int(cx.size), flush=True)
    fig_dir = ROOT / "figures"
    pdf = fig_dir / "cluster_38857678.pdf"
    save_paper_figure(fig, pdf, dpi=200)
    print("WROTE", pdf, flush=True)

if __name__ == "__main__":
    render(ROOT/"supplement/cluster_illustration.npz")
