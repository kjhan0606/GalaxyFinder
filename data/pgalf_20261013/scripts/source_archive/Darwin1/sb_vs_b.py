#!/usr/bin/env python3
"""Surface brightness at the stellar-FoF boundary, one halo at a time.

Star records are 120 bytes (NCHEM=3). The image uses every star in the
halo. Satellites (other galaxies with at least 100 stars) are masked.
Linking length scales from the b=0.2 particle mass. Ages invert the
RAMSES supercomoving time tau = int_1^a da/(a^3 E(a)), which matches
info time at snapshot 967.
"""
import os
os.environ.setdefault("SPS_HOME", "/home/kjhan/fsps_data")

import ctypes
import struct
from pathlib import Path

import numpy as np
from numba import njit
from scipy.spatial import cKDTree
from astropy.cosmology import FlatLambdaCDM
import astropy.units as u

HI = struct.Struct("<6i11d")
SI = struct.Struct("<5i4x11d")
DM, GAS, SINK, STAR = 120, 112, 168, 120
OM = 0.315
H0 = 67.4
ZSUN_SIM = 0.0134
COSMO = FlatLambdaCDM(H0=H0, Om0=OM)
# 1/H0 in Gyr. 977.792221 Gyr per (km/s/Mpc)^-1.
THUB = 977.792221 / H0

ROOT = Path("/gpfs/kjhan/Others/Jaehyun/2026.Apr")


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


def age_gyr(tp, a_snap, a_grid, tau_grid):
    """Star age at the snapshot. tp is RAMSES tau at birth."""
    a_b = np.interp(tp, tau_grid, a_grid)
    # proper-time integral tabulated
    # age = THUB * int_{a_b}^{a_snap} da/(a E(a))
    # precompute cumulative proper time from 0.02
    return a_b


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


def payload(n):
    return int(n[0]) * DM + int(n[1]) * GAS + int(n[2]) * SINK + int(n[3]) * STAR


def read_stars(f, counts):
    f.seek(int(counts[0]) * DM + int(counts[1]) * GAS + int(counts[2]) * SINK, 1)
    n = int(counts[3])
    if n == 0:
        return None
    blob = f.read(n * STAR)
    if len(blob) != n * STAR:
        raise RuntimeError("short star read")
    raw = np.frombuffer(blob, dtype=np.uint8).reshape(n, STAR)
    def f8(a, b):
        return raw[:, a:b].copy().view("<f8").reshape(n)
    def i4(a, b):
        return raw[:, a:b].copy().view("<i4").reshape(n)
    return {
        "x": f8(0, 8), "y": f8(8, 16), "z": f8(16, 24),
        "mass": f8(48, 56), "id": i4(56, 60),
        "tp": f8(72, 80), "zp": f8(80, 88), "mass0": f8(112, 120),
    }


def skip_payload(f, counts):
    f.seek(payload(counts), 1)


def iter_halos(path, background=False):
    """Yield (nstar, mstar, pieces) where pieces is a list of (npstar, center, stars or None)."""
    with open(path, "rb") as f:
        while True:
            if not background:
                buf = f.read(HI.size)
                if not buf:
                    return
                if len(buf) != HI.size:
                    raise RuntimeError("truncated halo")
                nsub, ndm, nstar, nsink, ngas, npall, totm, mdm, mgas, msink, mstar, x, y, z, vx, vy, vz = HI.unpack(buf)
            else:
                nsub = 1
                nstar = None
                mstar = None
            pieces = []
            for _ in range(nsub):
                rec = f.read(SI.size)
                if len(rec) != SI.size:
                    raise RuntimeError("truncated sub")
                npdm, npgas, npsink, npstar, npall_s, totm_s, mdm_s, mgas_s, msink_s, mstar_s, xs, ys, zs, vxs, vys, vzs = SI.unpack(rec)
                stars = read_stars(f, (npdm, npgas, npsink, npstar))
                pieces.append((npstar, (xs, ys, zs), stars, mstar_s))
            yield nstar, mstar, pieces


@njit
def _find(parent, i):
    while parent[i] != i:
        parent[i] = parent[parent[i]]
        i = parent[i]
    return i


@njit
def _unite_cells(x, y, z, link, ix, iy, iz, uniq, start, count, order, nx, ny):
    """Friend-of-friends using only occupied cells. uniq is sorted cell keys."""
    n = x.shape[0]
    parent = np.empty(n, np.int64)
    for i in range(n):
        parent[i] = i
    nu = uniq.shape[0]
    for i in range(n):
        for dz in range(-1, 2):
            az = iz[i] + dz
            for dy in range(-1, 2):
                ay = iy[i] + dy
                for dx in range(-1, 2):
                    ax = ix[i] + dx
                    key = np.int64(ax) + np.int64(nx) * (np.int64(ay) + np.int64(ny) * np.int64(az))
                    c = np.searchsorted(uniq, key)
                    if c >= nu or uniq[c] != key:
                        continue
                    jpos = start[c]
                    jend = jpos + count[c]
                    while jpos < jend:
                        j = order[jpos]
                        jpos += 1
                        if j <= i:
                            continue
                        r2 = (x[i] - x[j]) ** 2 + (y[i] - y[j]) ** 2 + (z[i] - z[j]) ** 2
                        lim = 0.5 * (link[i] + link[j])
                        if r2 < lim * lim:
                            ri = _find(parent, i)
                            rj = _find(parent, j)
                            if ri != rj:
                                if ri < rj:
                                    parent[rj] = ri
                                else:
                                    parent[ri] = rj
    root = np.empty(n, np.int64)
    for i in range(n):
        root[i] = _find(parent, i)
    return root


_OPFOF_SO = Path(__file__).resolve().parents[1] / "opFoF/variant/libopfof_sb.so"
_OPFOF = ctypes.CDLL(str(_OPFOF_SO))
_F64 = np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags="C_CONTIGUOUS")
_U64 = np.ctypeslib.ndpointer(dtype=np.uint64, ndim=1, flags="C_CONTIGUOUS")
_OPFOF.opfof_sb_labels.argtypes = [
    ctypes.c_size_t, _F64, _F64, _F64, _F64, _U64,
]
_OPFOF.opfof_sb_labels.restype = ctypes.c_int


def fof_largest(x, y, z, link):
    """Component label of each star from the opFoF octree walker."""
    n = int(x.shape[0])
    if n == 0:
        return np.zeros(0, np.int64)
    link = np.ascontiguousarray(link, dtype=np.float64)
    if float(np.max(link)) <= 0.0:
        return np.arange(n, dtype=np.int64)
    x = np.ascontiguousarray(x, dtype=np.float64)
    y = np.ascontiguousarray(y, dtype=np.float64)
    z = np.ascontiguousarray(z, dtype=np.float64)
    labels = np.empty(n, dtype=np.uint64)
    status = _OPFOF.opfof_sb_labels(n, x, y, z, link, labels)
    if status:
        raise RuntimeError("opFoF adapter failed with status %s for n=%s" % (status, n))
    return labels.astype(np.int64)


def link02(mass, b):
    return b * (mass / 2.7755e11 / OM) ** (1.0 / 3.0)


def scan_offsets(path, background=False):
    rows = []
    with open(path, "rb") as f:
        ih = 0
        while True:
            start = f.tell()
            if not background:
                buf = f.read(HI.size)
                if not buf:
                    break
                if len(buf) != HI.size:
                    raise RuntimeError("truncated halo")
                nsub, ndm, nstar = struct.unpack_from("<3i", buf)
                mstar = struct.unpack_from("<11d", buf, 24)[4]
            else:
                nsub, nstar, mstar = 1, None, None
                if f.read(1) == b"":
                    break
                f.seek(-1, 1)
            for _ in range(nsub):
                rec = f.read(SI.size)
                if len(rec) != SI.size:
                    raise RuntimeError("truncated sub at halo %d" % ih)
                npdm, npgas, npsink, npstar = struct.unpack_from("<4i", rec)
                f.seek(payload((npdm, npgas, npsink, npstar)), 1)
            rows.append((ih, nstar, mstar, start))
            ih += 1
            if ih % 100000 == 0:
                print("  offsets", ih, flush=True)
    return rows


def load_halo(path, offset, background=False):
    with open(path, "rb") as f:
        f.seek(offset)
        if not background:
            buf = f.read(HI.size)
            nsub = struct.unpack_from("<i", buf)[0]
        else:
            nsub = 1
        pieces = []
        for _ in range(nsub):
            rec = f.read(SI.size)
            npdm, npgas, npsink, npstar, npall, totm, mdm, mgas, msink, mstar, xs, ys, zs, vx, vy, vz = SI.unpack(rec)
            stars = read_stars(f, (npdm, npgas, npsink, npstar))
            pieces.append((npstar, np.array([xs, ys, zs]), stars, mstar))
    return pieces


def brightest(pieces):
    """Return the galaxy with most star particles. Legacy name, not luminosity."""
    k = int(np.argmax([p[0] for p in pieces]))
    return pieces[k]


def heaviest_mask(root, mass):
    roots, inv = np.unique(root, return_inverse=True)
    tot = np.bincount(inv, weights=mass)
    return root == roots[int(np.argmax(tot))]


def nucleus_index(x, y, z):
    """Star with the smallest radius enclosing itself and 31 neighbours.

    The mass-weighted centre can sit between clumps. The finder keeps the
    FoF group grown from the density peak. cKDTree k=32 includes the
    queried particle itself among the returned neighbours.
    """
    n = int(x.shape[0])
    if n <= 1:
        return 0
    k = 32 if n > 32 else n
    pts = np.column_stack((np.ascontiguousarray(x), np.ascontiguousarray(y), np.ascontiguousarray(z)))
    dist, _ = cKDTree(pts).query(pts, k=k, workers=1)
    return int(np.argmin(dist[:, -1]))


def peak_mask(root, x, y, z, i_peak=None):
    if i_peak is None:
        i_peak = nucleus_index(x, y, z)
    return root == root[int(i_peak)]


def ssp_flux(z, cache):
    """Flux of 1 Msun formed, observer frame, shape (nage, nz, 2) for V and g."""
    import fsps
    cache.parent.mkdir(parents=True, exist_ok=True)
    if cache.exists():
        d = np.load(cache)
        return d["age"], d["logz"], d["flux"]
    ZSUN_FSPS = 0.0142
    age = np.geomspace(0.01, 13.5, 32)
    logz = np.linspace(np.log10(1e-4 / ZSUN_FSPS), np.log10(0.05 / ZSUN_FSPS), 12)
    sp = fsps.StellarPopulation(zcontinuous=1, imf_type=1, sfh=0, dust_type=0,
                                dust2=0.0, add_neb_emission=False)
    flux = np.empty((age.size, logz.size, 2))
    for j, lz in enumerate(logz):
        sp.params["logzsol"] = float(lz)
        for i, tage in enumerate(age):
            mag = sp.get_mags(tage=float(tage), redshift=float(z), bands=["v", "sdss_g"])
            flux[i, j] = 10.0 ** (-0.4 * mag)
        print("ssp Z", j + 1, "/", logz.size, flush=True)
    np.savez(cache, age=age, logz=logz, flux=flux, z=z)
    return age, logz, flux


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


def profile(R_pkpc, flux, edges):
    """Sum flux in annuli. flux is (n, 2). Returns mu or nan."""
    idx = np.digitize(R_pkpc, edges) - 1
    nbin = edges.size - 1
    mu = np.full((nbin, 2), np.nan)
    area_pkpc = np.pi * (edges[1:] ** 2 - edges[:-1] ** 2)
    for b in range(2):
        sums = np.bincount(idx[(idx >= 0) & (idx < nbin)], weights=flux[(idx >= 0) & (idx < nbin), b], minlength=nbin)
        # convert pkpc^2 to arcsec^2 later by the caller via a scale
        mu[:, b] = sums
    return mu, area_pkpc


def crossing(R, mu, thresh):
    """Outer isophote: first faint bin after the profile has been brighter."""
    seen_bright = False
    for i in range(mu.size):
        if not np.isfinite(mu[i]):
            continue
        if mu[i] <= thresh:
            seen_bright = True
        elif seen_bright:
            return float(R[i])
    return np.nan


def interp_mu(r, Rc, mu):
    ok = np.isfinite(mu) & np.isfinite(Rc) & (Rc > 0)
    if ok.sum() < 2 or not np.isfinite(r) or r <= 0:
        return np.nan
    rr = Rc[ok]
    yy = mu[ok]
    if r < rr[0] or r > rr[-1]:
        return np.nan
    return float(np.interp(np.log(r), np.log(rr), yy))


def solve_b(bs, mus, thresh):
    """b where mu crosses the isophote while getting fainter as b grows.

    The outermost such crossing is the boundary of the main body. A crossing
    from faint to bright is a clump inside the profile and is left out.
    """
    bs = np.asarray(bs, float)
    mus = np.asarray(mus, float)
    ok = np.isfinite(bs) & np.isfinite(mus)
    bs, mus = bs[ok], mus[ok]
    if bs.size < 2:
        return np.nan
    order = np.argsort(bs)
    bs, mus = bs[order], mus[order]
    cross = np.nan
    for i in range(bs.size - 1):
        y0, y1 = mus[i], mus[i + 1]
        if y1 == y0:
            continue
        if y0 <= thresh < y1:
            frac = (thresh - y0) / (y1 - y0)
            cross = float(bs[i] + frac * (bs[i + 1] - bs[i]))
    return cross


def measure(pieces, bg, a_snap, z, grids):
    br = brightest(pieces)
    cx, cy, cz = br[1]
    box = 44.171267
    h = H0 / 100.0
    chunks = []
    for p in pieces:
        st = p[2]
        if st is None:
            continue
        sat = (p is not br) and p[0] >= 100
        chunks.append((st, np.full(st["mass"].size, sat), np.full(st["mass"].size, p is br)))
    if bg is not None and bg[2] is not None:
        st = bg[2]
        chunks.append((st, np.zeros(st["mass"].size, bool), np.zeros(st["mass"].size, bool)))
    cat = {k: np.concatenate([c[0][k] for c in chunks]) for k in ("x", "y", "z", "mass", "mass0", "tp", "zp")}
    sat = np.concatenate([c[1] for c in chunks])
    is_br = np.concatenate([c[2] for c in chunks])
    def unwrap(q, c):
        d = q - c
        return d - box * np.round(d / box)
    dx = unwrap(cat["x"], cx)
    dy = unwrap(cat["y"], cy)
    dz = unwrap(cat["z"], cz)
    R = np.hypot(dx, dy) / h * a_snap * 1000.0  # proper kpc
    a_grid, tau_grid, t_cum, age_g, logz_g, flux_g = grids
    # mass0 is Msun/h. FSPS magnitudes are for 1 solar mass formed.
    h_little = H0 / 100.0
    m0 = np.clip(cat["mass0"], 0, None) / h_little
    good = m0 > 0
    flux = np.zeros((cat["x"].size, 2))
    flux[good] = star_flux(cat["tp"][good], cat["zp"][good], m0[good], a_snap,
                           a_grid, tau_grid, t_cum, age_g, logz_g, flux_g)
    age_br = proper_age(np.clip(np.interp(cat["tp"][is_br], tau_grid, a_grid), 0.02, a_snap),
                        a_snap, a_grid, t_cum)
    mratio = cat["mass"][is_br] / np.clip(cat["mass0"][is_br], 1e-30, None)
    print("  age", float(age_br.min()), float(np.median(age_br)), float(age_br.max()),
          "zp", float(np.median(cat["zp"][is_br])),
          "m/m0", float(np.median(mratio)), flush=True)
    edges = np.geomspace(0.05, 3000.0, 64)
    Rc = np.sqrt(edges[:-1] * edges[1:])
    use = ~sat
    summed, area_pk = profile(R[use], flux[use], edges)
    DM = COSMO.comoving_transverse_distance(z).to_value(u.Mpc)
    # 1 proper kpc subtends arcsec: proper_kpc / (a * D_M_kpc) * 206265, D_M is comoving
    # theta(") = R_proper_Mpc / D_A_Mpc * 206265, D_A = a * D_M
    as_per_pkpc = (1e-3) / (a_snap * DM) * 206265.0
    area_as = area_pk * as_per_pkpc ** 2
    mu = np.full_like(summed, np.nan)
    for b in range(2):
        pos = summed[:, b] > 0
        mu[pos, b] = -2.5 * np.log10(summed[pos, b] / area_as[pos])
    R26 = crossing(Rc, mu[:, 0], 26.5)
    R27 = crossing(Rc, mu[:, 1], 27.0)
    # FoF only the b=0.2 brightest members
    xb, yb, zb = dx[is_br], dy[is_br], dz[is_br]
    mb = cat["mass"][is_br]
    fb = flux[is_br]
    Rb = R[is_br]
    i_peak = nucleus_index(xb, yb, zb)
    rows = []
    for bpar in (0.02, 0.03, 0.05, 0.08, 0.10, 0.12, 0.15, 0.18, 0.20):
        link = link02(mb, bpar)
        root = fof_largest(xb, yb, zb, link)
        keep = peak_mask(root, xb, yb, zb, i_peak)
        if keep.sum() < 10:
            rows.append((bpar, np.nan, np.nan, np.nan, np.nan, int(keep.sum())))
            continue
        r90s = []
        mus = []
        for band in (0, 1):
            light = fb[keep, band]
            rr = Rb[keep]
            if light.sum() <= 0:
                r90s.append(np.nan)
                mus.append(np.nan)
                continue
            order = np.argsort(rr)
            cume = np.cumsum(light[order])
            r90 = float(rr[order][np.searchsorted(cume, 0.9 * cume[-1])])
            r90s.append(r90)
            mus.append(interp_mu(r90, Rc, mu[:, band]))
        rows.append((bpar, r90s[0], mus[0], r90s[1], mus[1], int(keep.sum())))
    return rows, R26, R27, int(is_br.sum()), int(cat["x"].size)


B_GRID = (0.02, 0.03, 0.05, 0.08, 0.10, 0.12, 0.15, 0.18, 0.20)
# Rest-frame 26.5/27 isophotes are larger. The overlap was still rising at 0.20.
# 0.50 is a few times the mean density for equal-mass stars, and a peak there
# stays censored. The observer-frame grid is unchanged.
REST_B_GRID = B_GRID + (0.25, 0.30, 0.40, 0.50)
OUT = ROOT / "Darwin1/analysis_sb"
# Fixed mock image. These do not change with b or with the line of sight.
PIX_AS = 0.2
STAMP_AS = 180.0
KERNEL_AS = 80.0
FWHM_AS = 1.0
MOFFAT_BETA = 4.765
# DECam DR9–DR10 outer aureole, nmgy arcsec^-2 per nmgy of star, outside 7–8".
# V is not a DECam band; the r-band amplitude is the proxy.
WING_A = (0.00033, 0.00045)
THRESHOLDS = (26.5, 27.0)
# Rest-frame zero point. One pc^2 at 10 pc covers (206265/10)^2 arcsec^2.
# mu = -2.5 log10(flux per pc^2) + ABS_SB_ZERO, with flux = 10^(-0.4 M_abs).
ABS_SB_ZERO = float(2.5 * np.log10((206265.0 / 10.0) ** 2))
BAND_NAME = ("V", "g")
# Line of sight, then the two sky axes in (x, y, z).
LOS_AXES = ((1, 2), (2, 0), (0, 1))
LOS_NAME = ("x", "y", "z")
# Same stellar-mass window as the light end of the z=0.227 top 100.
MLO, MHI = 5.0e10, 1.20e11
# step, a, GALFIND, background, and optional nstar list (None = mass window)
STEPS = {
    "00967": (
        0.814714654255096,
        ROOT / "Darwin1/top100/1/GALFIND.DATA.00967",
        ROOT / "Darwin1/top100/1/background_ptl.00967",
        (986605, 1235506, 1501981, 1630920, 1902782, 2059707),
    ),
    "00081": (
        0.24987299740314484,
        ROOT / "Darwin1_recv/s00081/FoF_Data/FoF.00081/GALFIND.DATA.00081",
        ROOT / "Darwin1_recv/s00081/FoF_Data/FoF.00081/background_ptl.00081",
        None,
    ),
    "00178": (
        0.3344480097293854,
        ROOT / "Darwin1_recv/s00178/FoF_Data/FoF.00178/GALFIND.DATA.00178",
        ROOT / "Darwin1_recv/s00178/FoF_Data/FoF.00178/background_ptl.00178",
        None,
    ),
    "00475": (
        0.49948298931121826,
        ROOT / "Darwin1_recv/s00475/FoF_Data/FoF.00475/GALFIND.DATA.00475",
        ROOT / "Darwin1_recv/s00475/FoF_Data/FoF.00475/background_ptl.00475",
        None,
    ),
    "00786": (
        0.6665250062942505,
        ROOT / "Darwin1_recv/s00786/FoF_Data/FoF.00786/GALFIND.DATA.00786",
        ROOT / "Darwin1_recv/s00786/FoF_Data/FoF.00786/background_ptl.00786",
        None,
    ),
    "01001": (
        0.8400269746780396,
        ROOT / "Darwin1_recv/s01001/FoF_Data/FoF.01001/GALFIND.DATA.01001",
        ROOT / "Darwin1_recv/s01001/FoF_Data/FoF.01001/background_ptl.01001",
        None,
    ),
}


def cached_offsets(path, background, tag):
    cache = OUT / f"offsets_{tag}{'_bg' if background else ''}.npz"
    if cache.exists():
        d = np.load(cache)
        return list(zip(d["ih"].tolist(), d["nstar"].tolist(), d["mstar"].tolist(), d["off"].tolist()))
    print("scan", path, flush=True)
    rows = scan_offsets(path, background=background)
    ih = np.array([r[0] for r in rows], np.int64)
    nstar = np.array([(-1 if r[1] is None else r[1]) for r in rows], np.int64)
    mstar = np.array([(np.nan if r[2] is None else r[2]) for r in rows], np.float64)
    off = np.array([r[3] for r in rows], np.int64)
    cache.parent.mkdir(parents=True, exist_ok=True)
    np.savez(cache, ih=ih, nstar=nstar, mstar=mstar, off=off)
    print("cached", cache, "n", len(rows), flush=True)
    return rows


# Group mass (HaloInfo totm, Msun/h). At z<=0.23 a cluster is >=3e13.
# Earlier, a fixed 3e13 cut is almost empty, so the progenitors are >=1e13.
CLUSTER_MIN = {
    "00967": 3.0e13,
    "01001": 3.0e13,
    "00786": 1.0e13,
    "00475": 1.0e13,
    "00178": 1.0e13,
    "00081": 1.0e13,
}


def choose_cluster(step, gpath, rows):
    """Halos above the cluster group-mass cut, smallest star count first."""
    thr = CLUSTER_MIN[step]
    chosen = []
    if step == "00967":
        mass_by_n = {}
        text = (ROOT / "Darwin1/top100/top100_halos.txt").read_text().splitlines()
        for ln in text[1:]:
            if not ln.strip():
                continue
            col = ln.split()
            mass_by_n[int(col[3])] = float(col[2])
        for r in rows:
            totm = mass_by_n.get(int(r[1]))
            if totm is not None and totm >= thr:
                chosen.append((r, totm))
    else:
        with open(gpath, "rb") as f:
            for r in rows:
                if r[2] is None or not (r[2] >= 3.0e10):
                    continue
                f.seek(int(r[3]))
                buf = f.read(32)
                if len(buf) < 32:
                    continue
                totm = struct.unpack_from("<d", buf, 24)[0]
                if totm >= thr:
                    chosen.append((r, totm))
    chosen.sort(key=lambda item: item[0][1])
    return chosen, thr


def choose_halos(rows, want):
    if want is not None:
        by = {r[1]: r for r in rows}
        return [by[n] for n in want if n in by]
    cand = [r for r in rows if r[1] is not None and 5e5 <= r[1] <= 2.5e6 and MLO <= r[2] <= MHI]
    if len(cand) < 4:
        cand = [r for r in rows if r[1] is not None and 5e5 <= r[1] <= 3.0e6 and 3.0e10 <= r[2] <= 2.0e11]
    cand = sorted(cand, key=lambda r: r[2])
    if len(cand) <= 6:
        return cand
    pick = np.linspace(0, len(cand) - 1, 6).astype(int)
    return [cand[int(i)] for i in pick]


def _selftest_fof():
    x = np.array([0.0, 0.001, 0.1])
    y = np.zeros(3)
    z = np.zeros(3)
    link = np.array([0.002, 0.002, 0.002])
    root = fof_largest(x, y, z, link)
    assert root[0] == root[1] and root[0] != root[2]
    link = np.array([0.01, 0.001, 0.001])
    x = np.array([0.0, 0.005, 0.05])
    root = fof_largest(x, y, z, link)
    assert root[0] == root[1] and root[0] != root[2]


def check_fof():
    """Occupied-cell FoF at b=0.10 must reproduce the finder catalogue."""
    _selftest_fof()
    a_snap = 0.814714654255096
    a_grid, tau_grid, t_cum = build_time(a_snap)
    tau = float(np.interp(a_snap, a_grid, tau_grid))
    print("tau(a_snap)", tau, flush=True)
    if abs(tau + 0.23877) > 2e-3:
        raise SystemExit("tau(a_snap) is not -0.2388")
    g02 = ROOT / "Darwin1/top100/1/GALFIND.DATA.00967"
    g10 = ROOT / "Darwin1/top100/b0.10/GALFIND.DATA.00967"
    off02, off10 = 74541785392, 71033149728
    print("load catalogues", flush=True)
    p02 = load_halo(g02, off02)
    p10 = load_halo(g10, off10)
    br = brightest(p02)
    br10 = brightest(p10)
    st, st10 = br[2], br10[2]
    print("nstar b0.2", st["mass"].size, "catalogue b0.10", st10["mass"].size, flush=True)
    box = 44.171267
    cx, cy, cz = br[1]
    def unwrap(q, c):
        d = q - c
        return d - box * np.round(d / box)
    x, y, z = unwrap(st["x"], cx), unwrap(st["y"], cy), unwrap(st["z"], cz)
    link = link02(st["mass"], 0.10)
    root = fof_largest(x, y, z, link)
    cat = np.sort(st10["id"])
    ok = False
    for name, mask in (("peak", peak_mask(root, x, y, z)), ("heavy", heaviest_mask(root, st["mass"]))):
        got = np.sort(st["id"][mask])
        same = got.size == cat.size and np.array_equal(got, cat)
        print(name, "n", int(got.size), "equal", bool(same), flush=True)
        ok = ok or same
    if not ok:
        raise SystemExit("FoF does not reproduce the b=0.10 catalogue")
    print("CHECK_DONE", flush=True)


def run_step(step):
    a_snap, gpath, bpath, want = STEPS[step]
    z = 1.0 / a_snap - 1.0
    print("step", step, "z", z, flush=True)
    a_grid, tau_grid, t_cum = build_time(a_snap)
    print("tau(a_snap)", float(np.interp(a_snap, a_grid, tau_grid)), flush=True)
    cache = OUT / f"fsps_z{z:.3f}.npz"
    age_g, logz_g, flux_g = ssp_flux(z, cache)
    grids = (a_grid, tau_grid, t_cum, age_g, logz_g, flux_g)
    halos = cached_offsets(gpath, False, step)
    backs = cached_offsets(bpath, True, step)
    chosen = choose_halos(halos, want)
    print("chosen", [(r[0], r[1], r[2]) for r in chosen], flush=True)
    outp = OUT / f"sb_vs_b_{step}.txt"
    done = set()
    lines = []
    if outp.exists():
        lines = outp.read_text().splitlines()
        for line in lines:
            if line.startswith("H "):
                done.add(int(line.split()[2]))
    else:
        lines = [
            "# observer-frame SB at R90 of the FoF component that contains the stellar density peak",
            "# image is every star in the halo; satellites with >=100 stars are masked; area is kept",
            "# light = (mass0/h) * FSPS(age, Z=zp), Zsun_FSPS=0.0142, bands v and sdss_g, no dust/PSF/sky",
            "# H z nstar mstar b_V26.5 b_g27 R26.5 R27 n_bright",
            "# R z nstar b r90_V mu_V r90_g mu_g n_group",
        ]
    for ih, nstar, mstar, off in chosen:
        if nstar in done:
            print("skip", nstar, flush=True)
            continue
        print("halo", ih, nstar, mstar, flush=True)
        pieces = load_halo(gpath, off)
        bg = None
        if ih < len(backs):
            bg = load_halo(bpath, backs[ih][3], background=True)[0]
        rows, R26, R27, nbright, nall = measure(pieces, bg, a_snap, z, grids)
        bV = solve_b([r[0] for r in rows], [r[2] for r in rows], 26.5)
        bg_ = solve_b([r[0] for r in rows], [r[4] for r in rows], 27.0)
        print(" R26.5", R26, "R27", R27, "bV", bV, "bg", bg_, "nbright", nbright, "nall", nall, flush=True)
        lines.append(f"H {z:.4f} {nstar} {mstar:.6e} {bV:.4f} {bg_:.4f} {R26:.3f} {R27:.3f} {nbright}")
        for bpar, rV, muV, rG, mug, ng in rows:
            line = f"R {z:.4f} {nstar} {bpar:.2f} {rV:.3f} {muV:.3f} {rG:.3f} {mug:.3f} {ng}"
            lines.append(line)
            print(" ", line, flush=True)
        outp.write_text("\n".join(lines) + "\n")
        del pieces, bg
    print("wrote", outp, flush=True)


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


def _half_width(row, pix):
    """FWHM from a peak-centred 1D cut. row[0] is the peak."""
    half = row[0] / 2.0
    below = np.flatnonzero(row < half)
    if below.size == 0:
        return np.nan
    k = int(below[0])
    if k == 0:
        return np.nan
    y0, y1 = row[k - 1], row[k]
    frac = 0.0 if y1 == y0 else (half - y0) / (y1 - y0)
    return 2.0 * (k - 1 + frac) * pix


def check_psf():
    kern, core_sum, raw_sum = psf_kernel(PIX_AS, FWHM_AS, MOFFAT_BETA, WING_A[1], KERNEL_AS)
    fwhm = _half_width(kern[kern.shape[0] // 2, kern.shape[1] // 2:], PIX_AS)
    print("psf core", core_sum, "with_wing", raw_sum, "fwhm", fwhm, "n", kern.shape[0], flush=True)
    if abs(kern.sum() - 1.0) > 1e-8:
        raise SystemExit("PSF kernel is not normalized")
    if not (0.90 <= fwhm <= 1.10):
        raise SystemExit("PSF FWHM is not 1 arcsec")
    if not (raw_sum > core_sum and 0.002 < (raw_sum - core_sum) < 0.05):
        raise SystemExit("PSF wing integral is unexpected")
    print("PSF_OK", flush=True)
    return kern


def parent_stars(pieces, bg, sat_min=100):
    """Stars that enter the image and the linking.

    Satellites are the other catalogue galaxies with at least sat_min stars.
    On a stellar-FoF-off catalogue those are the dark-matter satellites.
    Background stars stay in. The centre belongs to the galaxy with the
    largest star-particle count, without a luminosity ranking.
    """
    br = brightest(pieces)
    centre = (float(br[1][0]), float(br[1][1]), float(br[1][2]))
    chunks = []
    n_sat = 0
    n_obj = 0
    for p in pieces:
        st = p[2]
        if st is None:
            continue
        if p is not br and int(p[0]) >= sat_min:
            n_sat += int(st["mass"].size)
            n_obj += 1
            continue
        chunks.append(st)
    if bg is not None and bg[2] is not None:
        chunks.append(bg[2])
    if not chunks:
        raise RuntimeError("halo has no stars after the satellite mask")
    keys = ("x", "y", "z", "mass", "mass0", "tp", "zp")
    cat = {k: np.concatenate([c[k] for c in chunks]) for k in keys}
    return cat, centre, n_sat, n_obj


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


def best_on_grid(bs, score):
    score = np.asarray(score, float)
    ok = np.isfinite(score)
    if not np.any(ok):
        return np.nan, np.nan, 1
    i = int(np.argmax(np.where(ok, score, -np.inf)))
    edge = int(i == 0 or i == bs.size - 1)
    return float(bs[i]), float(score[i]), edge


def measure_contour(pieces, bg, a_snap, z, grids, rest=False, b_grid=None):
    """Match each b to the fixed isophote on three axes.

    The surface-brightness map is every non-satellite star, convolved once.
    b is applied here, to that same star list. The score is the light-weighted
    overlap of the density-peak group with the region inside the isophote.
    rest keeps the angular pixels and the PSF, and thresholds rest-frame
    flux per proper pc^2. The flux grid must then be the redshift-0 cache.
    """
    from scipy.signal import fftconvolve
    cat, (cx, cy, cz), n_sat, n_obj = parent_stars(pieces, bg)
    box = 44.171267
    h = H0 / 100.0

    def unwrap(q, c):
        d = q - c
        return d - box * np.round(d / box)

    dx = unwrap(cat["x"], cx)
    dy = unwrap(cat["y"], cy)
    dz = unwrap(cat["z"], cz)
    mass = cat["mass"]
    a_grid, tau_grid, t_cum, age_g, logz_g, flux_g = grids
    m0 = np.clip(cat["mass0"], 0, None) / h
    good = m0 > 0
    flux = np.zeros((dx.size, 2))
    flux[good] = star_flux(cat["tp"][good], cat["zp"][good], m0[good], a_snap,
                           a_grid, tau_grid, t_cum, age_g, logz_g, flux_g)
    age = proper_age(np.clip(np.interp(cat["tp"], tau_grid, a_grid), 0.02, a_snap),
                     a_snap, a_grid, t_cum)
    print("  age", float(np.min(age)), float(np.median(age)), float(np.max(age)),
          "n_parent", int(dx.size), "n_sat", n_sat, "n_sat_obj", n_obj, flush=True)
    del cat, age
    DM = COSMO.comoving_transverse_distance(z).to_value(u.Mpc)
    as_per_pkpc = 1e-3 / (a_snap * DM) * 206265.0
    scale = (a_snap * 1000.0 / h) * as_per_pkpc
    print("  arcsec/pkpc", float(as_per_pkpc), "stamp_pkpc", STAMP_AS / as_per_pkpc, flush=True)
    rest_area = None
    if rest:
        pc_side = PIX_AS * 1000.0 / as_per_pkpc
        rest_area = pc_side * pc_side
        print("  rest-frame pc/pix", float(pc_side), "mu_zero", ABS_SB_ZERO, flush=True)
    n = int(round(2.0 * STAMP_AS / PIX_AS))
    kernels = []
    for amplitude in WING_A:
        kern, _, _ = psf_kernel(PIX_AS, FWHM_AS, MOFFAT_BETA, amplitude, KERNEL_AS)
        kernels.append(kern)
    coords = (dx, dy, dz)
    inside = np.zeros((3, 2, dx.size), np.bool_)
    areas = np.zeros((3, 2))
    for ia, (j, k) in enumerate(LOS_AXES):
        ix = np.rint(coords[j] * scale / PIX_AS).astype(np.int32) + n // 2
        iy = np.rint(coords[k] * scale / PIX_AS).astype(np.int32) + n // 2
        on = (ix >= 0) & (ix < n) & (iy >= 0) & (iy < n)
        for ib in range(2):
            img = bin_flux(ix, iy, flux[:, ib], n)
            img = fftconvolve(img, kernels[ib], mode="same")
            img[img < 0.0] = 0.0
            if rest:
                region, area = central_isophote(
                    img, THRESHOLDS[ib], PIX_AS, ABS_SB_ZERO, rest_area)
            else:
                region, area = central_isophote(img, THRESHOLDS[ib], PIX_AS)
            areas[ia, ib] = area
            inside[ia, ib, on] = region[iy[on], ix[on]]
            req = np.sqrt(area / np.pi) if area > 0 else 0.0
            print("  iso", LOS_NAME[ia], BAND_NAME[ib], "area", area, "r_eq", req, flush=True)
            del img, region
    i_peak = nucleus_index(dx, dy, dz)
    if b_grid is None:
        b_grid = REST_B_GRID if rest else B_GRID
    rows = []
    for bpar in b_grid:
        link = link02(mass, bpar)
        root = fof_largest(dx, dy, dz, link)
        keep = peak_mask(root, dx, dy, dz, i_peak)
        ng = int(keep.sum())
        print("  b", bpar, "n", ng, flush=True)
        for ia in range(3):
            for ib in range(2):
                w = flux[:, ib]
                inn = inside[ia, ib]
                both = keep & inn
                lg = float(w[keep].sum())
                li = float(w[inn].sum())
                lb = float(w[both].sum())
                den = lg + li - lb
                iou = lb / den if den > 0.0 else np.nan
                rec = lb / li if li > 0.0 else np.nan
                pur = lb / lg if lg > 0.0 else np.nan
                rows.append((bpar, ia, ib, iou, rec, pur, ng, areas[ia, ib]))
    return rows, int(dx.size), n_sat


def run_contour(step, cluster=False, rest=False, only_nstar=None):
    a_snap, gpath, bpath, want = STEPS[step]
    tag = step
    if step == "00967":
        # Stellar linking was skipped. Halo membership is the parent set.
        gpath = ROOT / "Darwin1/top100/0/GALFIND.DATA.00967"
        bpath = ROOT / "Darwin1/top100/0/background_ptl.00967"
        tag = "00967_off"
    z = 1.0 / a_snap - 1.0
    if rest:
        mode = "rest"
    elif cluster:
        mode = "cluster"
    else:
        mode = "contour"
    print(mode, step, "z", z, gpath, flush=True)
    check_psf()
    a_grid, tau_grid, t_cum = build_time(a_snap)
    if rest:
        cache = OUT / "fsps_rest_z0.npz"
        age_g, logz_g, flux_g = ssp_flux(0.0, cache)
    else:
        cache = OUT / f"fsps_z{z:.3f}.npz"
        age_g, logz_g, flux_g = ssp_flux(z, cache)
    grids = (a_grid, tau_grid, t_cum, age_g, logz_g, flux_g)
    halos = cached_offsets(gpath, False, tag)
    backs = cached_offsets(bpath, True, tag)
    if cluster:
        picked, thr = choose_cluster(step, gpath, halos)
        chosen = [item[0] for item in picked]
        totm_of = {int(item[0][1]): item[1] for item in picked}
        print("cluster cut", thr, "n", len(chosen), flush=True)
    else:
        chosen = choose_halos(halos, want)
        totm_of = {}
    print("chosen", [(r[0], r[1], r[2]) for r in chosen], flush=True)
    if only_nstar is not None:
        chosen = [r for r in chosen if int(r[1]) == int(only_nstar)]
        print("only nstar", int(only_nstar), "n", len(chosen), flush=True)
    if not chosen:
        raise SystemExit("no halos matched in %s" % gpath)
    if rest:
        outp = OUT / f"sb_rest_{step}.txt"
    elif cluster:
        outp = OUT / f"sb_cluster_{step}.txt"
    else:
        outp = OUT / f"sb_contour_{step}.txt"
    catalogue = outp
    if only_nstar is not None:
        part_dir = OUT / "parts"
        part_dir.mkdir(parents=True, exist_ok=True)
        outp = part_dir / ("%s_%d.txt" % (catalogue.stem, int(only_nstar)))
    b_grid = REST_B_GRID if rest else B_GRID
    done = set()
    lines = []
    resume_from = catalogue if only_nstar is not None else outp
    if resume_from.exists():
        lines = resume_from.read_text().splitlines()
        top = float(b_grid[-1])
        have_top = set()
        present = set()
        for line in lines:
            if line.startswith("I "):
                parts = line.split()
                if abs(float(parts[5]) - top) < 1e-6:
                    have_top.add(int(parts[2]))
            elif line.startswith("S "):
                present.add(int(line.split()[2]))
        redo = present - have_top
        if redo:
            kept = []
            for line in lines:
                if line[:2] in ("P ", "A ", "I ", "C ", "S "):
                    if int(line.split()[2]) in redo:
                        continue
                kept.append(line)
            lines = kept
        done = have_top
        if only_nstar is not None:
            if int(only_nstar) in done:
                print("skip", int(only_nstar), flush=True)
                print("wrote", catalogue, flush=True)
                return
            lines = ["# fragment %s nstar %d" % (catalogue.name, int(only_nstar))]
    elif rest:
        lines = [
            "# rest-frame surface density. Same angular pixels and the same PSF as sb_cluster.",
            "# FSPS redshift=0 absolute magnitudes. Flux of 1 Msun formed is 10^(-0.4 M). mass0 is divided by h.",
            "# pixel proper area (pc^2) = (0.2 arcsec * 1000 / arcsec_per_proper_kpc)^2",
            "# mu = -2.5 log10(flux_per_pc2) + 2.5 log10((206265/10)^2). Thresholds 26.5 V and 27 g.",
            "# IoU weights are this rest-frame flux. Observer-frame FSPS caches are not read.",
            "# A area is the isophote area in arcsec^2. The threshold uses proper pc^2.",
            "# PSF Moffat FWHM=1 arcsec beta=4.765 + DECam wing A/r^2, A_V=0.00033 (r proxy) A_g=0.00045",
            "# wing ramps 7-8 arcsec, kernel to 80 arcsec, renormalized to 1; pixel 0.2 arcsec; no sky",
            "# satellites: other catalogue galaxies with >=100 stars; background stays",
            "# cluster sample: z=0.227 and z=0.19 have group mass >= 3e13 Msun/h;"
            " earlier redshifts >= 1e13 Msun/h",
            "# P z nstar mstar n_parent n_sat group_mass",
            "# A z nstar los band area_arcsec2",
            "# I z nstar los band b iou recall purity n_group",
            "# C z nstar los band b_best iou edge",
            "# S z nstar band b_med iou_med n_edge",
        ]
    else:
        lines = [
            "# fixed all-star image, PSF once, b applied at analysis, three lines of sight",
            "# z=0.227 uses stellar-FoF-off top100/0; other steps reassemble the b=0.2 halo stars",
            "# satellites: other catalogue galaxies with >=100 stars; background stays",
            "# PSF Moffat FWHM=1 arcsec beta=4.765 + DECam wing A/r^2, A_V=0.00033 (r proxy) A_g=0.00045",
            "# wing ramps 7-8 arcsec, kernel to 80 arcsec, renormalized to 1; pixel 0.2 arcsec; no sky",
            "# IoU = L(group and inside) / (L(group)+L(inside)-L(both)); light is each star's flux",
            "# inside = central component of the map brighter than 26.5 V or 27 sdss_g",
            "# axes x,y,z are the line of sight; they are not stacked",
            "# cluster sample: z=0.227 and z=0.19 have group mass >= 3e13 Msun/h;"
            " earlier redshifts >= 1e13 Msun/h",
            "# P z nstar mstar n_parent n_sat group_mass",
            "# A z nstar los band area_arcsec2",
            "# I z nstar los band b iou recall purity n_group",
            "# C z nstar los band b_best iou edge",
            "# S z nstar band b_med iou_med n_edge",
        ]
    for ih, nstar, mstar, off in chosen:
        if nstar in done:
            print("skip", nstar, flush=True)
            continue
        totm = totm_of.get(int(nstar), np.nan)
        print("halo", ih, nstar, mstar, "M", totm, flush=True)
        pieces = load_halo(gpath, off)
        bg = None
        if ih < len(backs):
            bg = load_halo(bpath, backs[ih][3], background=True)[0]
        rows, n_parent, n_sat = measure_contour(
            pieces, bg, a_snap, z, grids, rest=rest, b_grid=b_grid)
        del pieces, bg
        lines.append(f"P {z:.4f} {nstar} {mstar:.6e} {n_parent} {n_sat} {totm:.6e}")
        seen_area = set()
        by = {}
        for bpar, ia, ib, iou, rec, pur, ng, area in rows:
            by.setdefault((ia, ib), []).append((bpar, iou, rec, pur, ng))
            if (ia, ib) not in seen_area:
                seen_area.add((ia, ib))
                lines.append(f"A {z:.4f} {nstar} {LOS_NAME[ia]} {BAND_NAME[ib]} {area:.3f}")
        bmed = []
        for (ia, ib), items in sorted(by.items()):
            bs = np.array([t[0] for t in items])
            ious = np.array([t[1] for t in items])
            b_best, iou_best, edge = best_on_grid(bs, ious)
            bmed.append((ib, b_best, iou_best, edge))
            lines.append(
                f"C {z:.4f} {nstar} {LOS_NAME[ia]} {BAND_NAME[ib]} {b_best:.4f} {iou_best:.4f} {edge}"
            )
            for bpar, iou, rec, pur, ng in items:
                line = (f"I {z:.4f} {nstar} {LOS_NAME[ia]} {BAND_NAME[ib]} {bpar:.2f} "
                        f"{iou:.4f} {rec:.4f} {pur:.4f} {ng}")
                lines.append(line)
                print(" ", line, flush=True)
        for ib, name in enumerate(BAND_NAME):
            vals = [t for t in bmed if t[0] == ib]
            bb = np.array([t[1] for t in vals], float)
            ii = np.array([t[2] for t in vals], float)
            n_edge = int(np.sum([t[3] for t in vals]))
            line = (f"S {z:.4f} {nstar} {name} {np.nanmedian(bb):.4f} "
                    f"{np.nanmedian(ii):.4f} {n_edge}")
            lines.append(line)
            print(" ", line, flush=True)
        outp.write_text("\n".join(lines) + "\n")
    print("wrote", outp, flush=True)


def main():
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("steps", nargs="*", default=["00967"])
    ap.add_argument("--check", action="store_true")
    ap.add_argument("--psf", action="store_true")
    ap.add_argument("--contour", action="store_true")
    ap.add_argument("--cluster", action="store_true")
    ap.add_argument("--rest", action="store_true",
                    help="rest-frame surface density on the cluster sample")
    ap.add_argument("--nstar", type=int, default=None,
                    help="compute one halo and write a fragment under analysis_sb/parts")
    args = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    if args.psf:
        check_psf()
        return
    if args.check:
        check_fof()
        return
    if args.rest and (args.cluster or args.contour):
        raise SystemExit("--rest is separate from --cluster and --contour")
    if args.rest:
        for step in args.steps:
            run_contour(step, cluster=True, rest=True, only_nstar=args.nstar)
        return
    if args.cluster:
        for step in args.steps:
            run_contour(step, cluster=True, only_nstar=args.nstar)
        return
    if args.contour:
        for step in args.steps:
            run_contour(step, only_nstar=args.nstar)
        return
    for step in args.steps:
        run_step(step)


if __name__ == "__main__":
    main()
