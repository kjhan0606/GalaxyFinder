"""Numerical definitions extracted from the calculation sources.

Inputs are catalogue stellar arrays and the packaged SSP grids. The module
contains the image and score calculations, including the corrected boundary
restriction. Stellar FoF membership is supplied as a Boolean array.
"""
import numpy as np
from scipy import ndimage
from scipy.signal import fftconvolve
OM = .315
H0 = 67.4
THUB = 977.792221 / H0
STRUCT = np.ones((3, 3), np.uint8)

def E(a):
    return np.sqrt(OM / a**3 + (1.0 - OM))

def tau_table(a_max, n=4000):
    """tau(a) = int_1^a da'/(a'^3 E). Negative for a < 1."""
    a = np.linspace(0.02, a_max, n)
    grid = np.linspace(a_max, 1.0, 4000)
    tau_max = -np.trapezoid(1.0 / (grid**3 * E(grid)), grid)
    inv = 1.0 / (a**3 * E(a))
    integ = np.zeros(n)
    integ[1:] = np.cumsum(0.5 * (inv[:-1] + inv[1:]) * np.diff(a))
    tau = (tau_max - integ[-1]) + integ
    return a, tau

def proper_age(a_birth, a_snap, a_grid, t_cum):
    tb = np.interp(a_birth, a_grid, t_cum)
    ts = np.interp(a_snap, a_grid, t_cum)
    return np.clip(ts - tb, 1e-3, None)

def build_time(a_snap):
    a_grid, tau_grid = tau_table(a_snap)
    # cumulative proper time from a_grid[0]
    inv = 1.0 / (a_grid * E(a_grid))
    dt = np.diff(a_grid)
    mid = 0.5 * (inv[:-1] + inv[1:])
    t_cum = np.zeros_like(a_grid)
    t_cum[1:] = np.cumsum(mid * dt) * THUB
    # check tau(a_snap) later
    return a_grid, tau_grid, t_cum

def star_flux(tp, zp, mass0, a_snap, a_grid, tau_grid, t_cum, age_g, logz_g, flux_g):
    a_b = np.interp(tp, tau_grid, a_grid)
    age = proper_age(a_b, a_snap, a_grid, t_cum)
    ZSUN_FSPS = 0.0142
    logz = np.log10(np.clip(zp, 1e-5, 0.05) / ZSUN_FSPS)
    # bilinear in (log age, logz) of the log flux
    la = np.log(age_g)
    lstar = np.log(np.clip(age, age_g[0], age_g[-1]))
    x = np.interp(lstar, la, np.arange(la.size))
    y = np.interp(logz, logz_g, np.arange(logz_g.size))
    x0 = np.floor(x).astype(int)
    y0 = np.floor(y).astype(int)
    x0 = np.clip(x0, 0, la.size - 2)
    y0 = np.clip(y0, 0, logz_g.size - 2)
    tx = x - x0
    ty = y - y0
    out = np.empty((tp.size, 2))
    for b in range(2):
        f = np.log(flux_g[:, :, b])
        f00 = f[x0, y0]
        f10 = f[x0 + 1, y0]
        f01 = f[x0, y0 + 1]
        f11 = f[x0 + 1, y0 + 1]
        out[:, b] = mass0 * np.exp((1 - tx) * (1 - ty) * f00 + tx * (1 - ty) * f10
                                   + (1 - tx) * ty * f01 + tx * ty * f11)
    return out

def link02(mass, b):
    return b * (mass / 2.7755e11 / OM) ** (1.0 / 3.0)

def moffat_alpha(fwhm, beta):
    return fwhm / (2.0 * np.sqrt(2.0 ** (1.0 / beta) - 1.0))

def psf_kernel(pix, fwhm, beta, amplitude, rmax):
    """Moffat seeing plus a DECam A/r^2 wing, then renormalized to 1.

    The wing is zero inside 7 arcsec and rises linearly to full strength
    at 8 arcsec, matching the Legacy Surveys halo model. The Moffat is the
    analytic unit-integral profile sampled at pixel centres.
    """
    n = int(np.ceil(rmax / pix))
    yy, xx = np.mgrid[-n:n + 1, -n:n + 1]
    r = pix * np.hypot(xx, yy)
    alpha = moffat_alpha(fwhm, beta)
    shape = (beta - 1.0) / (np.pi * alpha * alpha) * (1.0 + (r / alpha) ** 2) ** (-beta)
    core = shape * (pix * pix)
    wing = np.zeros_like(core)
    outer = r >= 7.0
    ramp = np.clip((r - 7.0) / 1.0, 0.0, 1.0)
    wing[outer] = ramp[outer] * amplitude / (r[outer] ** 2) * (pix * pix)
    raw = core + wing
    total = float(raw.sum())
    return raw / total, float(core.sum()), total

def bin_flux(ix, iy, weight, n):
    img = np.zeros((n, n), np.float64)
    m = (ix >= 0) & (ix < n) & (iy >= 0) & (iy < n) & np.isfinite(weight) & (weight > 0)
    if np.any(m):
        np.add.at(img, (iy[m], ix[m]), weight[m])
    return img

def central_isophote(img, thresh, pix, mag_offset=0.0, area=None):
    """Pixels brighter than thresh that connect to the stamp centre.

    The default divisor is pix*pix, so the map is flux per arcsec^2.
    area replaces that divisor. mag_offset is added after -2.5 log10.
    The returned area is always the region size in arcsec^2.
    """
    from scipy import ndimage
    if area is None:
        per_as = img / (pix * pix)
        bright = per_as > 0
        mu = np.full(img.shape, np.inf)
        mu[bright] = -2.5 * np.log10(per_as[bright])
    else:
        per = img / area
        bright = per > 0
        mu = np.full(img.shape, np.inf)
        mu[bright] = -2.5 * np.log10(per[bright]) + mag_offset
    mask = mu <= thresh
    labels, _ = ndimage.label(mask, structure=np.ones((3, 3), np.uint8))
    cy = cx = img.shape[0] // 2
    if labels[cy, cx] == 0:
        ys, xs = np.nonzero(mask)
        if ys.size == 0:
            return np.zeros(img.shape, np.bool_), 0.0
        j = int(np.argmin((ys - cy) ** 2 + (xs - cx) ** 2))
        keep = labels[ys[j], xs[j]]
    else:
        keep = labels[cy, cx]
    region = labels == keep
    return region, float(region.sum()) * pix * pix

def pixels(coord_a, coord_b, scale, n, pix=0.2):
    ix = np.rint(coord_a * scale / pix).astype(np.int32) + n // 2
    iy = np.rint(coord_b * scale / pix).astype(np.int32) + n // 2
    on = (ix >= 0) & (ix < n) & (iy >= 0) & (iy < n)
    return ix, iy, on

def sample_inside(region, ix, iy, on):
    inside = np.zeros(ix.shape[0], np.bool_)
    inside[on] = region[iy[on], ix[on]]
    return inside

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
