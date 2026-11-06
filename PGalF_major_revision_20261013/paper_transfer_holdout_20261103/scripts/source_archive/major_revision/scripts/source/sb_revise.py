#!/usr/bin/env python3
"""Observer-frame revision measurements. Does not rewrite sb_cluster tables.

The baseline image, isophote, and friends-of-friends scores follow
sb_vs_b.measure_contour. A halo is kept only when those scores reproduce
the published I lines. Threshold and sky re-threshold the same convolved
image. The wing-off and satellite-in images are extra convolutions. Thinning
rebuilds the image and the friends-of-friends group.
"""
import sys
from pathlib import Path

import numpy as np
from scipy import ndimage
from scipy.signal import fftconvolve

ROOT = Path("/gpfs/kjhan/Others/Jaehyun/2026.Apr")
sys.path.insert(0, str(Path(__file__).resolve().parent))
import sb_vs_b as sb

REV = sb.OUT / "revision"
B_GRID = np.asarray(sb.B_GRID, float)
LOS = ("x", "y", "z")
BAND = ("V", "g")
# Representative failures, plus two interior successes with modest star counts.
FIG_NSTARS = {
    10033180, 7062464, 17494810, 5557744, 9003659, 4640948,
    27762214, 18705275, 11111575, 11019776, 4927969,
    4211038, 6770450,
}
STRUCT = np.ones((3, 3), np.uint8)


def catalogue_paths(step):
    a_snap, gpath, bpath, _want = sb.STEPS[step]
    tag = step
    if step == "00967":
        gpath = ROOT / "Darwin1/top100/0/GALFIND.DATA.00967"
        bpath = ROOT / "Darwin1/top100/0/background_ptl.00967"
        tag = "00967_off"
    z = 1.0 / a_snap - 1.0
    return a_snap, z, gpath, bpath, tag


def published_lines(step, nstar):
    found = []
    path = sb.OUT / f"sb_cluster_{step}.txt"
    for line in path.read_text().splitlines():
        if line[:2] in ("A ", "I "):
            if int(line.split()[2]) == int(nstar):
                found.append(line)
    return found


def stack(chunks):
    keys = ("x", "y", "z", "mass", "mass0", "tp", "zp")
    return {k: np.concatenate([c[k] for c in chunks]) for k in keys}


def satellite_chunks(pieces):
    br = sb.brightest(pieces)
    out = []
    for p in pieces:
        if p[2] is None:
            continue
        if p is not br and int(p[0]) >= 100:
            out.append(p[2])
    return out


def unwrap_flux(cat, centre, a_snap, grids):
    """Same coordinates and fluxes as measure_contour."""
    cx, cy, cz = centre
    box = 44.171267
    h = sb.H0 / 100.0

    def unwrap(q, c):
        d = q - c
        return d - box * np.round(d / box)

    dx = unwrap(cat["x"], cx)
    dy = unwrap(cat["y"], cy)
    dz = unwrap(cat["z"], cz)
    mass = np.ascontiguousarray(cat["mass"], dtype=np.float64)
    a_grid, tau_grid, t_cum, age_g, logz_g, flux_g = grids
    m0 = np.clip(cat["mass0"], 0, None) / h
    good = m0 > 0
    flux = np.zeros((dx.size, 2))
    if np.any(good):
        flux[good] = sb.star_flux(
            cat["tp"][good], cat["zp"][good], m0[good], a_snap,
            a_grid, tau_grid, t_cum, age_g, logz_g, flux_g)
    return dx, dy, dz, mass, flux


def angular_scale(a_snap, z):
    from astropy import units as u
    dm = sb.COSMO.comoving_transverse_distance(z).to_value(u.Mpc)
    as_per_pkpc = 1e-3 / (a_snap * dm) * 206265.0
    h = sb.H0 / 100.0
    scale = (a_snap * 1000.0 / h) * as_per_pkpc
    n = int(round(2.0 * sb.STAMP_AS / sb.PIX_AS))
    return float(scale), n


def pixels(coord_a, coord_b, scale, n):
    ix = np.rint(coord_a * scale / sb.PIX_AS).astype(np.int32) + n // 2
    iy = np.rint(coord_b * scale / sb.PIX_AS).astype(np.int32) + n // 2
    on = (ix >= 0) & (ix < n) & (iy >= 0) & (iy < n)
    return ix, iy, on


def sample_inside(region, ix, iy, on):
    inside = np.zeros(ix.shape[0], np.bool_)
    inside[on] = region[iy[on], ix[on]]
    return inside


def mu_labels(img, thresh):
    per = img / (sb.PIX_AS * sb.PIX_AS)
    bright = per > 0
    mu = np.full(img.shape, np.inf)
    mu[bright] = -2.5 * np.log10(per[bright])
    mask = mu <= thresh
    labels, nlab = ndimage.label(mask, structure=STRUCT)
    return mu, mask, labels, int(nlab)


def component_at(labels, mask, py, px):
    ny, nx = labels.shape
    off = 0
    if not (0 <= py < ny and 0 <= px < nx):
        off = 1
        py = int(np.clip(py, 0, ny - 1))
        px = int(np.clip(px, 0, nx - 1))
    if labels[py, px] != 0:
        keep = int(labels[py, px])
        fallback = 0
    else:
        ys, xs = np.nonzero(mask)
        if ys.size == 0:
            return 0, 1, off
        j = int(np.argmin((ys - py) ** 2 + (xs - px) ** 2))
        keep = int(labels[ys[j], xs[j]])
        fallback = 1
    return keep, fallback, off


def region_of(labels, keep, pix):
    if keep == 0:
        return np.zeros(labels.shape, np.bool_), 0.0
    region = labels == keep
    return region, float(region.sum()) * pix * pix


def light_scores(w, keep, inside):
    both = keep & inside
    lg = float(w[keep].sum())
    li = float(w[inside].sum())
    lb = float(w[both].sum())
    den = lg + li - lb
    iou = lb / den if den > 0.0 else np.nan
    rec = lb / li if li > 0.0 else np.nan
    pur = lb / lg if lg > 0.0 else np.nan
    return iou, rec, pur


def ring_of(region):
    if region.sum() == 0:
        return region
    eroded = ndimage.binary_erosion(region, structure=STRUCT)
    dilated = ndimage.binary_dilation(region, structure=STRUCT)
    return (region & ~eroded) | (dilated & ~region)


def boundary_scores(w, keep, inside, in_ring):
    """Compare group and original image selection only within the boundary band.

    Earlier revision tables incorrectly compared the full group to the band.
    The corrected intersection is G & I & R and the union is (G | I) & R.
    Archived rev_*.txt files retain the earlier calculation and must not be
    interpreted as the output of this corrected function.
    """
    return light_scores(w[in_ring], keep[in_ring], inside[in_ring])


def best_curve(ious):
    return sb.best_on_grid(B_GRID, np.asarray(ious, float))


def fof_keeps(dx, dy, dz, mass):
    i_peak = sb.nucleus_index(dx, dy, dz)
    keeps = []
    for bpar in B_GRID:
        link = sb.link02(mass, bpar)
        root = sb.fof_largest(dx, dy, dz, link)
        keeps.append(sb.peak_mask(root, dx, dy, dz, i_peak))
        print("  b", float(bpar), "n", int(keeps[-1].sum()), flush=True)
    return keeps, i_peak


def convolve_band(raw, kernel):
    img = fftconvolve(raw, kernel, mode="same")
    img[img < 0.0] = 0.0
    return img


def project_raw(dx, dy, dz, flux, scale, n):
    coords = (dx, dy, dz)
    raws = []
    pixs = []
    for j, k in sb.LOS_AXES:
        ix, iy, on = pixels(coords[j], coords[k], scale, n)
        raws.append([sb.bin_flux(ix, iy, flux[:, ib], n) for ib in range(2)])
        pixs.append((ix, iy, on))
    return raws, pixs


def curves_for(flux, keeps, inside):
    """IoU, recall, purity at each grid b. inside is per band, shared across b."""
    out = []
    for ib in range(2):
        w = flux[:, ib]
        row = []
        for keep in keeps:
            row.append(light_scores(w, keep, inside[ib]))
        out.append(row)
    return out


def emit_cur(fh, tag, z, nstar, areas, curves, rings):
    """rings[los][band][b_index] is the ring IoU."""
    for ia, los in enumerate(LOS):
        for ib, band in enumerate(BAND):
            ious = [curves[ia][ib][t][0] for t in range(len(B_GRID))]
            b_best, iou_best, edge = best_curve(ious)
            i_best = int(np.nanargmax(np.where(np.isfinite(ious), ious, -np.inf))) if np.any(np.isfinite(ious)) else 0
            def at(b):
                k = int(np.argmin(np.abs(B_GRID - b)))
                return rings[ia][ib][k]
            vals = " ".join(f"{v:.6f}" if np.isfinite(v) else "nan" for v in ious)
            rb = rings[ia][ib][i_best]
            fh.write(
                f"CUR {tag} {z:.4f} {nstar} {los} {band} {areas[ia][ib]:.3f} "
                f"{b_best:.4f} {iou_best:.6f} {edge} {vals} "
                f"{rb:.6f} {at(0.10):.6f} {at(0.20):.6f}\n"
            )


def ring_ious(flux, keeps, ix, iy, on, region):
    ious = []
    for ib in range(2):
        inside = sample_inside(region, ix, iy, on)
        in_ring = sample_inside(ring_of(region), ix, iy, on)
        ious.append([boundary_scores(flux[:, ib], keep, inside, in_ring)[0] for keep in keeps])
    return ious


def write_radial(fh, z, nstar, sky, flux, keeps, scale, n, areas, b_list):
    for ia, los in enumerate(LOS):
        ca, cb = sky[ia][0], sky[ia][1]
        r_as = np.hypot(ca, cb) * scale
        for ib, band in enumerate(BAND):
            area = areas[ia][ib]
            req = np.sqrt(area / np.pi) if area > 0 else np.nan
            if not np.isfinite(req) or req <= 0:
                continue
            rn = r_as / req
            w = flux[:, ib]
            bins = np.digitize(rn, [0.5, 1.0, 1.5])
            for bpar in b_list:
                k = int(np.argmin(np.abs(B_GRID - bpar)))
                keep = keeps[k]
                inside = sky[ia][3][ib]
                lin_all = float(w[inside].sum())
                lgrp_all = float(w[keep].sum())
                for ibin in range(4):
                    m = bins == ibin
                    sel_in = m & inside
                    sel_g = m & keep
                    both = sel_in & keep
                    lin = float(w[sel_in].sum())
                    lg = float(w[sel_g].sum())
                    lb = float(w[both].sum())
                    rec = lb / lin if lin > 0 else np.nan
                    pur = lb / lg if lg > 0 else np.nan
                    fin = lin / lin_all if lin_all > 0 else np.nan
                    fgrp = lg / lgrp_all if lgrp_all > 0 else np.nan
                    fh.write(
                        f"RAD {z:.4f} {nstar} {los} {band} {B_GRID[k]:.2f} {ibin} "
                        f"{rec:.6f} {pur:.6f} {fin:.6f} {fgrp:.6f} "
                        f"{lb:.6e} {lin:.6e} {lg:.6e}\n"
                    )


def write_sb(fh, z, nstar, mu_at, on, flux, keeps, inside, thresh, b_list):
    """mu_at is inf off the stamp. Bins are the faint side of the isophote."""
    names = ("lt22", "b22_24", "b24_26", "b26_thr", "thr_p1", "faint", "off")
    for ia, los in enumerate(LOS):
        for ib, band in enumerate(BAND):
            mu = np.full(flux.shape[0], np.inf)
            mu[on[ia]] = mu_at[ia][ib][on[ia]]
            thr = thresh[ib]
            edges_mu = mu
            # bin id
            bid = np.full(mu.shape, 6, np.int8)  # off / zero
            finite = np.isfinite(edges_mu)
            bid[finite & (edges_mu < 22)] = 0
            bid[finite & (edges_mu >= 22) & (edges_mu < 24)] = 1
            bid[finite & (edges_mu >= 24) & (edges_mu < 26)] = 2
            bid[finite & (edges_mu >= 26) & (edges_mu <= thr)] = 3
            bid[finite & (edges_mu > thr) & (edges_mu <= thr + 1)] = 4
            bid[finite & (edges_mu > thr + 1)] = 5
            w = flux[:, ib]
            for bpar in b_list:
                k = int(np.argmin(np.abs(B_GRID - bpar)))
                keep = keeps[k]
                inn = inside[ia][ib]
                lin_all = float(w[inn].sum())
                lgrp_all = float(w[keep].sum())
                for ibin, name in enumerate(names):
                    m = bid == ibin
                    sel_in = m & inn
                    sel_g = m & keep
                    both = sel_in & keep
                    lin = float(w[sel_in].sum())
                    lg = float(w[sel_g].sum())
                    lb = float(w[both].sum())
                    rec = lb / lin if lin > 0 else np.nan
                    pur = lb / lg if lg > 0 else np.nan
                    fin = lin / lin_all if lin_all > 0 else np.nan
                    fgrp = lg / lgrp_all if lgrp_all > 0 else np.nan
                    fh.write(
                        f"SB {z:.4f} {nstar} {los} {band} {B_GRID[k]:.2f} {name} "
                        f"{rec:.6f} {pur:.6f} {fin:.6f} {fgrp:.6f} "
                        f"{lb:.6e} {lin:.6e} {lg:.6e}\n"
                    )


def variant_pack(flux, keeps, pixs, images, thresh_pair, tag):
    """Re-threshold stored images. Returns areas, curves, ring ious per los."""
    areas = []
    insides = []
    ring_pack = []
    for ia in range(3):
        ix, iy, on = pixs[ia]
        area_b = []
        inside_b = []
        rings_b = []
        for ib in range(2):
            region, area = sb.central_isophote(images[ia][ib], thresh_pair[ib], sb.PIX_AS)
            inside_b.append(sample_inside(region, ix, iy, on))
            area_b.append(area)
            ring = ring_of(region)
            rings_b.append(ring)
        areas.append(area_b)
        insides.append(inside_b)
        # ring IoU per band per b, after both bands exist
        per_band = []
        for ib in range(2):
            inside_r = sample_inside(rings_b[ib], ix, iy, on)
            per_band.append([boundary_scores(flux[:, ib], keep, inside_b[ib], inside_r)[0] for keep in keeps])
        ring_pack.append(per_band)
    # curves: [los][band][b] -> scores, built from insides
    curves = []
    for ia in range(3):
        curves.append(curves_for(flux, keeps, insides[ia]))
    return areas, curves, ring_pack


def save_figure(path, img, region, peak_region, ix, iy, on, keep, ix_peak, iy_peak, n,
                title, zoom=None):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    per = img / (sb.PIX_AS * sb.PIX_AS)
    mu = np.full(img.shape, np.nan)
    bright = per > 0
    mu[bright] = -2.5 * np.log10(per[bright])
    half = 0.5 * n * sb.PIX_AS
    extent = [-half, half, -half, half]
    fig, axes = plt.subplots(1, 2, figsize=(9.2, 4.6), constrained_layout=True)
    for ax in axes:
        ax.imshow(mu, origin="lower", cmap="gray_r", vmin=20, vmax=28, extent=extent)
        if region.any():
            ax.contour(region.astype(float), levels=[0.5], colors="#D97757",
                       linewidths=0.7, extent=extent, origin="lower")
        if peak_region is not None and peak_region.any():
            ax.contour(peak_region.astype(float), levels=[0.5], colors="#4EC9B0",
                       linewidths=0.6, extent=extent, origin="lower")
        ax.plot(0, 0, marker="+", color="#3898EC", ms=8, mew=1.2)
        if 0 <= ix_peak < n and 0 <= iy_peak < n:
            ax.plot((ix_peak - n / 2) * sb.PIX_AS, (iy_peak - n / 2) * sb.PIX_AS,
                    marker="x", color="#FFC107", ms=7, mew=1.2)
        ax.set_aspect("equal")
    # Group members on the right panel only.
    m = keep & on
    if m.sum() > 6000:
        rng = np.random.default_rng(1)
        choose = np.flatnonzero(m)
        choose = rng.choice(choose, 6000, replace=False)
        m = np.zeros(keep.shape, np.bool_)
        m[choose] = True
    xs = (ix[m] - n / 2) * sb.PIX_AS
    ys = (iy[m] - n / 2) * sb.PIX_AS
    axes[1].scatter(xs, ys, s=1, c="#3898EC", alpha=0.35, linewidths=0, rasterized=True)
    axes[0].set_title("isophote")
    axes[1].set_title("density-peak group")
    if zoom is not None:
        for ax in axes:
            ax.set_xlim(-zoom, zoom)
            ax.set_ylim(-zoom, zoom)
    fig.suptitle(title, fontsize=9)
    fig.savefig(path, dpi=140)
    plt.close(fig)


def prepare_halo(step, nstar):
    a_snap, z, gpath, bpath, tag = catalogue_paths(step)
    print("revise", step, "z", z, "nstar", nstar, flush=True)
    sb.check_psf()
    a_grid, tau_grid, t_cum = sb.build_time(a_snap)
    cache = sb.OUT / f"fsps_z{z:.3f}.npz"
    age_g, logz_g, flux_g = sb.ssp_flux(z, cache)
    grids = (a_grid, tau_grid, t_cum, age_g, logz_g, flux_g)
    halos = sb.cached_offsets(gpath, False, tag)
    backs = sb.cached_offsets(bpath, True, tag)
    chosen = [r for r in halos if int(r[1]) == int(nstar)]
    if len(chosen) != 1:
        raise SystemExit("nstar %s matched %d halos in %s" % (nstar, len(chosen), gpath))
    ih, nstar, mstar, off = chosen[0]
    nstar = int(nstar)
    picked, _thr = sb.choose_cluster(step, gpath, [chosen[0]])
    totm = float(picked[0][1]) if picked else np.nan
    pieces = sb.load_halo(gpath, off)
    bg = None
    if ih < len(backs):
        bg = sb.load_halo(bpath, backs[ih][3], background=True)[0]
    cat, centre, n_sat, n_obj = sb.parent_stars(pieces, bg)
    sats = satellite_chunks(pieces)
    sat = stack(sats) if sats else None
    del pieces, bg, sats
    dx, dy, dz, mass, flux = unwrap_flux(cat, centre, a_snap, grids)
    del cat
    print("  n_parent", int(dx.size), "n_sat", n_sat, "n_sat_obj", n_obj,
          "M", totm, flush=True)
    mfin = mass[np.isfinite(mass) & (mass > 0)]
    ell = sb.link02(mfin, 0.10)
    print("  MASS med", float(np.median(mfin)), "ell10", float(np.median(ell)), flush=True)
    return {
        "a_snap": a_snap, "z": z, "nstar": nstar, "mstar": float(mstar),
        "totm": totm, "n_parent": int(dx.size), "n_sat": int(n_sat),
        "dx": dx, "dy": dy, "dz": dz, "mass": mass, "flux": flux,
        "sat": sat, "grids": grids, "centre": centre, "step": step,
        "m_med": float(np.median(mfin)), "ell_med": float(np.median(ell)),
        "m_p16": float(np.quantile(mfin, 0.16)), "m_p84": float(np.quantile(mfin, 0.84)),
        "ell_p16": float(np.quantile(ell, 0.16)), "ell_p84": float(np.quantile(ell, 0.84)),
    }


def run_one(step, nstar):
    REV.mkdir(parents=True, exist_ok=True)
    (REV / "figs").mkdir(parents=True, exist_ok=True)
    pub = published_lines(step, nstar)
    if not pub:
        raise SystemExit("no published lines for %s %s" % (step, nstar))
    H = prepare_halo(step, nstar)
    z = H["z"]
    nstar = H["nstar"]
    scale, n = angular_scale(H["a_snap"], z)
    print("  arcsec scale", scale, "n", n, flush=True)
    kernels = []
    for amplitude in sb.WING_A:
        kern, _, _ = sb.psf_kernel(sb.PIX_AS, sb.FWHM_AS, sb.MOFFAT_BETA, amplitude, sb.KERNEL_AS)
        kernels.append(kern)
    kern0, _, _ = sb.psf_kernel(sb.PIX_AS, sb.FWHM_AS, sb.MOFFAT_BETA, 0.0, sb.KERNEL_AS)
    raws, pixs = project_raw(H["dx"], H["dy"], H["dz"], H["flux"], scale, n)
    images = []
    for ia in range(3):
        images.append([convolve_band(raws[ia][ib], kernels[ib]) for ib in range(2)])
    keeps, i_peak = fof_keeps(H["dx"], H["dy"], H["dz"], H["mass"])
    # Baseline isophotes, with the centre / nearest-centre rule of the paper.
    base_areas = []
    base_inside = []
    base_region = []
    base_mu = []
    base_labels = []
    base_mask = []
    diag = []
    peak_region = []
    coords = (H["dx"], H["dy"], H["dz"])
    ixp = []
    iyp = []
    for ia, (j, k) in enumerate(sb.LOS_AXES):
        ix, iy, on = pixs[ia]
        area_b, inside_b, region_b, mu_b = [], [], [], []
        labels_b, mask_b, diag_b, preg_b = [], [], [], []
        px = int(np.rint(coords[j][i_peak] * scale / sb.PIX_AS)) + n // 2
        py = int(np.rint(coords[k][i_peak] * scale / sb.PIX_AS)) + n // 2
        ixp.append(px)
        iyp.append(py)
        for ib in range(2):
            mu, mask, labels, nlab = mu_labels(images[ia][ib], sb.THRESHOLDS[ib])
            cid, fallback, off = component_at(labels, mask, n // 2, n // 2)
            region, area = region_of(labels, cid, sb.PIX_AS)
            pid, pfall, poff = component_at(labels, mask, py, px)
            preg, parea = region_of(labels, pid, sb.PIX_AS)
            inside_b.append(sample_inside(region, ix, iy, on))
            area_b.append(area)
            region_b.append(region)
            mu_b.append(mu)
            labels_b.append(labels)
            mask_b.append(mask)
            same = int(cid == pid and cid != 0)
            sep = sb.PIX_AS * float(np.hypot(px - n / 2, py - n / 2))
            diag_b.append((nlab, area, parea, same, fallback, off, poff, sep, cid, pid))
            preg_b.append(preg)
            print("  iso", LOS[ia], BAND[ib], "area", area, "fallback", fallback,
                  "same_peak", same, "sep", sep, flush=True)
        base_areas.append(area_b)
        base_inside.append(inside_b)
        base_region.append(region_b)
        base_mu.append(mu_b)
        base_labels.append(labels_b)
        base_mask.append(mask_b)
        diag.append(diag_b)
        peak_region.append(preg_b)
    # Scores.
    base_curves = []
    for ia in range(3):
        base_curves.append(curves_for(H["flux"], keeps, base_inside[ia]))
    # Compare to the catalogue before any variant.
    got = []
    for ia, los in enumerate(LOS):
        for ib, band in enumerate(BAND):
            got.append(f"A {z:.4f} {nstar} {los} {band} {base_areas[ia][ib]:.3f}")
            for t, bpar in enumerate(B_GRID):
                iou, rec, pur = base_curves[ia][ib][t]
                ng = int(keeps[t].sum())
                got.append(
                    f"I {z:.4f} {nstar} {los} {band} {bpar:.2f} "
                    f"{iou:.4f} {rec:.4f} {pur:.4f} {ng}"
                )
    pub_set = set(pub)
    got_set = set(got)
    missing = sorted(pub_set - got_set)
    extra = sorted(got_set - pub_set)
    if missing or extra:
        print("BASELINE_FAIL", step, nstar, flush=True)
        for line in missing[:12]:
            print("  missing", line, flush=True)
        for line in extra[:12]:
            print("  extra", line, flush=True)
        raise SystemExit(2)
    print("BASELINE_MATCH", step, nstar, flush=True)

    outp = REV / f"rev_{step}_{nstar}.txt"
    with outp.open("w") as fh:
        fh.write("# revision measurements. Baseline reproduces sb_cluster I and A lines.\n")
        fh.write("# sky 27 is a 1 percent residual of a 22 mag/arcsec^2 sky; "
                 "sky 29 is 0.16 percent. Added after the kernel clip.\n")
        fh.write("# wing0 uses amplitude 0 and is renormalised to one.\n")
        fh.write("# thinning keeps each drawn star's mass and luminosity. "
                 "The link formula is unchanged. Nested uniforms, "
                 "seed sequence [20261007, nstar].\n")
        fh.write("# CUR tag z nstar los band area b_best iou_best edge "
                 "iou at the nine grid nodes, ring IoU at b_best, 0.10, 0.20\n")
        fh.write(
            f"MASS {z:.4f} {nstar} {H['n_parent']} {H['m_med']:.8e} "
            f"{H['m_p16']:.8e} {H['m_p84']:.8e} {H['ell_med']:.8e} "
            f"{H['ell_p16']:.8e} {H['ell_p84']:.8e} {H['totm']:.8e} {H['mstar']:.8e}\n"
        )
        for ia, los in enumerate(LOS):
            for ib, band in enumerate(BAND):
                nlab, area, parea, same, fallback, off, poff, sep, cid, pid = diag[ia][ib]
                fh.write(
                    f"CMP {z:.4f} {nstar} {los} {band} {nlab} {area:.3f} {parea:.3f} "
                    f"{same} {fallback} {off} {poff} {sep:.4f}\n"
                )
        # Positional versus photometric light, baseline mask.
        # Satellite light is convolved separately and is not in the baseline image.
        sat = H["sat"]
        sat_flux = None
        sat_pix = None
        sat_raw = None
        if sat is not None:
            sdx, sdy, sdz, smass, sat_flux = unwrap_flux(sat, H["centre"], H["a_snap"], H["grids"])
            del smass
            sat_raw, sat_pix = project_raw(sdx, sdy, sdz, sat_flux, scale, n)
            H["sat_coord"] = (sdx, sdy, sdz)
        for ia in range(3):
            for ib in range(2):
                photo_in = float(images[ia][ib][base_region[ia][ib]].sum())
                photo_tot = float(images[ia][ib].sum())
                pos_in = float(H["flux"][base_inside[ia][ib], ib].sum())
                pos_tot = float(H["flux"][:, ib].sum())
                if sat_flux is None:
                    spin = spot = sphin = spht = 0.0
                else:
                    six, siy, son = sat_pix[ia]
                    sin = sample_inside(base_region[ia][ib], six, siy, son)
                    spin = float(sat_flux[sin, ib].sum())
                    spot = float(sat_flux[:, ib].sum())
                    simg = convolve_band(sat_raw[ia][ib], kernels[ib])
                    sphin = float(simg[base_region[ia][ib]].sum())
                    spht = float(simg.sum())
                    sat_raw[ia][ib] = simg  # store the convolved satellite image
                fh.write(
                    f"PHO {z:.4f} {nstar} {LOS[ia]} {BAND[ib]} "
                    f"{pos_in:.8e} {pos_tot:.8e} {photo_in:.8e} {photo_tot:.8e} "
                    f"{spin:.8e} {spot:.8e} {sphin:.8e} {spht:.8e}\n"
                )
        # Baseline curve and rings.
        base_rings = []
        for ia in range(3):
            ix, iy, on = pixs[ia]
            per_band = []
            for ib in range(2):
                inside_r = sample_inside(ring_of(base_region[ia][ib]), ix, iy, on)
                per_band.append([boundary_scores(H["flux"][:, ib], keep, base_inside[ia][ib], inside_r)[0] for keep in keeps])
            base_rings.append(per_band)
        emit_cur(fh, "base", z, nstar, base_areas, base_curves, base_rings)
        # Radius and surface-brightness bins at the optimum, at 0.10, and at 0.20.
        sky = []
        for ia, (j, k) in enumerate(sb.LOS_AXES):
            sky.append((coords[j], coords[k], None, base_inside[ia]))
        b_wanted = []
        for ia in range(3):
            for ib in range(2):
                ious = [base_curves[ia][ib][t][0] for t in range(len(B_GRID))]
                b_best, _, _ = best_curve(ious)
                for bpar in (b_best, 0.10, 0.20):
                    if np.isfinite(bpar):
                        b_wanted.append(float(bpar))
        b_list = sorted(set(np.round(b_wanted, 5)))
        write_radial(fh, z, nstar, sky, H["flux"], keeps, scale, n, base_areas, b_list)
        mu_at = []
        on_list = []
        for ia in range(3):
            ix, iy, on = pixs[ia]
            mu_at.append([
                _mu_on_stars(base_mu[ia][ib], ix, iy, on, H["flux"].shape[0])
                for ib in range(2)
            ])
            on_list.append(on)
        # write_sb expects mu sampled on stars. Rebuild simply here.
        _write_sb_stars(fh, z, nstar, pixs, base_mu, H["flux"], keeps, base_inside,
                        sb.THRESHOLDS, b_list, n)
        # Figures for the listed halos, at the g and V images.
        if nstar in FIG_NSTARS:
            for ia, los in enumerate(LOS):
                for ib, band in enumerate(BAND):
                    ious = [base_curves[ia][ib][t][0] for t in range(len(B_GRID))]
                    b_best, iou_best, edge = best_curve(ious)
                    k = int(np.nanargmax(np.where(np.isfinite(ious), ious, -np.inf)))
                    title = (f"z={z:.3f}  nstar={nstar}  LOS {los}  {band}  "
                             f"area={base_areas[ia][ib]:.2f}  b={b_best:.2f}  IoU={iou_best:.3f}")
                    figp = REV / "figs" / f"{nstar}_{los}_{band}.png"
                    try:
                        save_figure(
                            figp, images[ia][ib], base_region[ia][ib], peak_region[ia][ib],
                            pixs[ia][0], pixs[ia][1], pixs[ia][2], keeps[k],
                            ixp[ia], iyp[ia], n, title)
                        area = base_areas[ia][ib]
                        if area < 5.0:
                            zoom = 5.0 if area <= 0 else max(3.0, 8.0 * np.sqrt(area / np.pi))
                            save_figure(
                                REV / "figs" / f"{nstar}_{los}_{band}_zoom.png",
                                images[ia][ib], base_region[ia][ib], peak_region[ia][ib],
                                pixs[ia][0], pixs[ia][1], pixs[ia][2], keeps[k],
                                ixp[ia], iyp[ia], n, title, zoom=zoom)
                        print("  fig", figp, flush=True)
                    except Exception as exc:
                        print("  FIG_FAIL", exc, flush=True)
        # Thresholds. Only the moved band is a new measurement; the other band
        # is re-thresholded at its baseline value and will match base.
        thresh_sets = [
            ("thV26.0", (26.0, sb.THRESHOLDS[1])),
            ("thV27.0", (27.0, sb.THRESHOLDS[1])),
            ("thg26.5", (sb.THRESHOLDS[0], 26.5)),
            ("thg27.5", (sb.THRESHOLDS[0], 27.5)),
        ]
        for tag, pair in thresh_sets:
            areas, curves, rings = variant_pack(H["flux"], keeps, pixs, images, pair, tag)
            emit_cur(fh, tag, z, nstar, areas, curves, rings)
            print("  variant", tag, flush=True)
        # Sky residual added to the clipped image, then the baseline threshold.
        for tag, mu_sky in (("sky29", 29.0), ("sky27", 27.0)):
            add = (sb.PIX_AS ** 2) * 10 ** (-0.4 * mu_sky)
            sky_images = [[img + add for img in los_imgs] for los_imgs in images]
            areas, curves, rings = variant_pack(
                H["flux"], keeps, pixs, sky_images, sb.THRESHOLDS, tag)
            emit_cur(fh, tag, z, nstar, areas, curves, rings)
            del sky_images
            print("  variant", tag, flush=True)
        # Wing removed.
        wing_images = []
        for ia in range(3):
            wing_images.append([convolve_band(raws[ia][ib], kern0) for ib in range(2)])
        areas, curves, rings = variant_pack(
            H["flux"], keeps, pixs, wing_images, sb.THRESHOLDS, "wing0")
        emit_cur(fh, "wing0", z, nstar, areas, curves, rings)
        print("  variant wing0", flush=True)
        del wing_images, raws
        # Component that contains the projected density peak.
        peak_areas = []
        peak_curves_inside = []
        for ia in range(3):
            ix, iy, on = pixs[ia]
            area_b, inside_b = [], []
            for ib in range(2):
                _nlab, _area, _parea, _same, _fb, _off, _poff, _sep, _cid, pid = diag[ia][ib]
                region, area = region_of(base_labels[ia][ib], pid, sb.PIX_AS)
                inside_b.append(sample_inside(region, ix, iy, on))
                area_b.append(area)
            peak_areas.append(area_b)
            peak_curves_inside.append(inside_b)
        peak_curves = [curves_for(H["flux"], keeps, peak_curves_inside[ia]) for ia in range(3)]
        peak_rings = []
        for ia in range(3):
            ix, iy, on = pixs[ia]
            per_band = []
            for ib in range(2):
                _nlab, _area, _parea, _same, _fb, _off, _poff, _sep, _cid, pid = diag[ia][ib]
                region, _area = region_of(base_labels[ia][ib], pid, sb.PIX_AS)
                inside_r = sample_inside(ring_of(region), ix, iy, on)
                per_band.append([boundary_scores(H["flux"][:, ib], keep, peak_curves_inside[ia][ib], inside_r)[0] for keep in keeps])
            peak_rings.append(per_band)
        emit_cur(fh, "peak", z, nstar, peak_areas, peak_curves, peak_rings)
        print("  variant peak", flush=True)
        # No nearest-centre fallback: an empty centre is an unresolved isophote.
        nf_areas, nf_inside = [], []
        for ia in range(3):
            ix, iy, on = pixs[ia]
            area_b, inside_b = [], []
            for ib in range(2):
                fallback = diag[ia][ib][4]
                if fallback:
                    region = np.zeros((n, n), np.bool_)
                    area = 0.0
                else:
                    region = base_region[ia][ib]
                    area = base_areas[ia][ib]
                inside_b.append(sample_inside(region, ix, iy, on))
                area_b.append(area)
            nf_areas.append(area_b)
            nf_inside.append(inside_b)
        nf_curves = [curves_for(H["flux"], keeps, nf_inside[ia]) for ia in range(3)]
        nf_rings = []
        for ia in range(3):
            ix, iy, on = pixs[ia]
            per_band = []
            for ib in range(2):
                if diag[ia][ib][4]:
                    per_band.append([np.nan] * len(B_GRID))
                else:
                    inside_r = sample_inside(ring_of(base_region[ia][ib]), ix, iy, on)
                    per_band.append([boundary_scores(H["flux"][:, ib], keep, nf_inside[ia][ib], inside_r)[0] for keep in keeps])
            nf_rings.append(per_band)
        emit_cur(fh, "nofb", z, nstar, nf_areas, nf_curves, nf_rings)
        print("  variant nofb", flush=True)
        # Satellites left in the image, still excluded from the link.
        if sat_raw is not None:
            both_images = []
            for ia in range(3):
                both_images.append([images[ia][ib] + sat_raw[ia][ib] for ib in range(2)])
            areas, curves, rings = variant_pack(
                H["flux"], keeps, pixs, both_images, sb.THRESHOLDS, "satin")
            emit_cur(fh, "satin", z, nstar, areas, curves, rings)
            # Photometric correspondence on the satellite-in isophote.
            for ia in range(3):
                for ib in range(2):
                    region, area = sb.central_isophote(
                        both_images[ia][ib], sb.THRESHOLDS[ib], sb.PIX_AS)
                    photo_in = float(both_images[ia][ib][region].sum())
                    photo_tot = float(both_images[ia][ib].sum())
                    ix, iy, on = pixs[ia]
                    inn = sample_inside(region, ix, iy, on)
                    pos_in = float(H["flux"][inn, ib].sum())
                    pos_tot = float(H["flux"][:, ib].sum())
                    six, siy, son = sat_pix[ia]
                    sin = sample_inside(region, six, siy, son)
                    spin = float(sat_flux[sin, ib].sum())
                    spot = float(sat_flux[:, ib].sum())
                    sphin = float(sat_raw[ia][ib][region].sum())
                    spht = float(sat_raw[ia][ib].sum())
                    fh.write(
                        f"PHU {z:.4f} {nstar} {LOS[ia]} {BAND[ib]} {area:.3f} "
                        f"{pos_in:.8e} {pos_tot:.8e} {photo_in:.8e} {photo_tot:.8e} "
                        f"{spin:.8e} {spot:.8e} {sphin:.8e} {spht:.8e}\n"
                    )
            del both_images
            print("  variant satin", flush=True)
        del images, base_region, base_mu, base_labels
        # Thinning. Nested keep, masses and luminosities unchanged.
        rng = np.random.default_rng(np.random.SeedSequence([20261007, nstar]))
        draw = rng.random(H["dx"].size)
        for tag, frac in (("thin2", 0.5), ("thin4", 0.25), ("thin8", 0.125)):
            m = draw < frac
            print("  thin", tag, "n", int(m.sum()), flush=True)
            tdx, tdy, tdz = H["dx"][m], H["dy"][m], H["dz"][m]
            tmass = H["mass"][m]
            tflux = H["flux"][m]
            traw, tpix = project_raw(tdx, tdy, tdz, tflux, scale, n)
            timages = []
            for ia in range(3):
                timages.append([convolve_band(traw[ia][ib], kernels[ib]) for ib in range(2)])
            tkeeps, _tp = fof_keeps(tdx, tdy, tdz, tmass)
            areas, curves, rings = variant_pack(
                tflux, tkeeps, tpix, timages, sb.THRESHOLDS, tag)
            emit_cur(fh, tag, z, nstar, areas, curves, rings)
            mfin = tmass[np.isfinite(tmass) & (tmass > 0)]
            ell = sb.link02(mfin, 0.10)
            fh.write(
                f"THIN {z:.4f} {nstar} {tag} {int(m.sum())} "
                f"{float(np.median(mfin)):.8e} {float(np.median(ell)):.8e}\n"
            )
            del traw, timages, tkeeps, tflux
            print("  variant", tag, flush=True)
        fh.write(f"DONE {z:.4f} {nstar}\n")
    print("WROTE", outp, flush=True)


def _mu_on_stars(mu_img, ix, iy, on, n):
    out = np.full(n, np.inf)
    out[on] = mu_img[iy[on], ix[on]]
    return out


def _region_from_area(img, thresh):
    region, _area = sb.central_isophote(img, thresh, sb.PIX_AS)
    return region


def _write_sb_stars(fh, z, nstar, pixs, mu_img, flux, keeps, inside, thresh, b_list, n):
    names = ("lt22", "b22_24", "b24_26", "b26_thr", "thr_p1", "faint", "off")
    for ia, los in enumerate(LOS):
        ix, iy, on = pixs[ia]
        for ib, band in enumerate(BAND):
            mu = np.full(flux.shape[0], np.inf)
            mu[on] = mu_img[ia][ib][iy[on], ix[on]]
            thr = thresh[ib]
            bid = np.full(mu.shape, 6, np.int8)
            finite = np.isfinite(mu)
            bid[finite & (mu < 22)] = 0
            bid[finite & (mu >= 22) & (mu < 24)] = 1
            bid[finite & (mu >= 24) & (mu < 26)] = 2
            bid[finite & (mu >= 26) & (mu <= thr)] = 3
            bid[finite & (mu > thr) & (mu <= thr + 1)] = 4
            bid[finite & (mu > thr + 1)] = 5
            w = flux[:, ib]
            for bpar in b_list:
                k = int(np.argmin(np.abs(B_GRID - bpar)))
                keep = keeps[k]
                inn = inside[ia][ib]
                lin_all = float(w[inn].sum())
                lgrp_all = float(w[keep].sum())
                for ibin, name in enumerate(names):
                    m = bid == ibin
                    sel_in = m & inn
                    sel_g = m & keep
                    both = sel_in & keep
                    lin = float(w[sel_in].sum())
                    lg = float(w[sel_g].sum())
                    lb = float(w[both].sum())
                    rec = lb / lin if lin > 0 else np.nan
                    pur = lb / lg if lg > 0 else np.nan
                    fin = lin / lin_all if lin_all > 0 else np.nan
                    fgrp = lg / lgrp_all if lgrp_all > 0 else np.nan
                    fh.write(
                        f"SB {z:.4f} {nstar} {los} {band} {B_GRID[k]:.2f} {name} "
                        f"{rec:.6f} {pur:.6f} {fin:.6f} {fgrp:.6f} "
                        f"{lb:.6e} {lin:.6e} {lg:.6e}\n"
                    )


def main():
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--step", required=True)
    ap.add_argument("--nstar", type=int, required=True)
    args = ap.parse_args()
    run_one(args.step, args.nstar)


if __name__ == "__main__":
    main()
