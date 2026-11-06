#!/usr/bin/env python3
"""Three checks that do not replace the adopted positional linking length.

The fine grid inserts 0.06, 0.07, 0.09, 0.11, 0.13, and 0.14 between the
published nodes and still ends at 0.02 and 0.20. The PSF score gives each
star the fraction of its own kernel that lands inside the adopted isophote.
Thinning keeps the published random subset, leaves each mass and luminosity
unchanged, and scores that subset against the original isophote. It does
not rebuild the image. Does not rewrite rev_*.txt or dist_*.txt.
"""
import sys
from pathlib import Path

import numpy as np
from scipy.signal import fftconvolve

REV = Path("/gpfs/kjhan/Others/Jaehyun/2026.Apr/Darwin1/analysis_sb/revision")
sys.path.insert(0, str(REV))
import sb_revise as rv
import sb_vs_b as sb

ORIG_B = np.asarray(rv.B_GRID, float)
FINE_B = np.asarray(
    (0.02, 0.03, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10, 0.11, 0.12, 0.13, 0.14, 0.15, 0.18, 0.20),
    float,
)
THIN = (("thin2", 0.5), ("thin4", 0.25), ("thin8", 0.125))


def fail(msg):
    print("FINE_FAIL", msg, flush=True)
    raise SystemExit(2)


def fmt(x):
    if not np.isfinite(x):
        return "nan"
    return "%.6f" % float(x)


def scores_frac(w, keep, frac):
    lg = float(w[keep].sum()) if np.any(keep) else 0.0
    li = float(np.sum(w * frac))
    lb = float(np.sum(w[keep] * frac[keep])) if np.any(keep) else 0.0
    den = lg + li - lb
    iou = lb / den if den > 0.0 else np.nan
    rec = lb / li if li > 0.0 else np.nan
    pur = lb / lg if lg > 0.0 else np.nan
    return iou, rec, pur


def self_test():
    w = np.array([1.0, 2.0, 3.0, 4.0])
    keep = np.array([True, True, False, True])
    inside = np.array([True, False, True, True])
    got = scores_frac(w, keep, inside.astype(float))
    ref = rv.light_scores(w, keep, inside)
    if any(abs(a - b) > 1e-12 for a, b in zip(got, ref)):
        fail("continuous score disagrees with the positional score %s %s" % (got, ref))
    kern = np.zeros((5, 5))
    kern[2, 2] = 0.4
    kern[2, 1] = kern[2, 3] = 0.2
    kern[1, 2] = kern[3, 2] = 0.1
    raw = np.zeros((9, 9))
    raw[4, 4] = 10.0
    img = fftconvolve(raw, kern, mode="same")
    region = np.zeros((9, 9), dtype=bool)
    region[4, 4] = True
    region[4, 5] = True
    frac = fftconvolve(region.astype(np.float64), kern, mode="same")
    photo = float(img[region].sum())
    pred = 10.0 * float(frac[4, 4])
    if abs(photo - pred) > 1e-9:
        fail("kernel fraction does not reproduce the aperture sum %s %s" % (photo, pred))
    if list(np.round(np.sort(np.concatenate([ORIG_B, [0.06, 0.07, 0.09, 0.11, 0.13, 0.14]])), 5)) != list(np.round(np.sort(FINE_B), 5)):
        fail("fine grid is not the published grid plus the six inserted nodes")
    print("SELF_TEST_OK", flush=True)


def load_base(step, nstar):
    path = REV / f"rev_{step}_{nstar}.txt"
    if not path.is_file():
        fail("missing %s" % path.name)
    text = path.read_text().splitlines()
    if not any(line.startswith("DONE ") for line in text):
        fail("rev file has no DONE line %s" % path.name)
    rows = {}
    for line in text:
        p = line.split()
        if len(p) < 19 or p[0] != "CUR" or p[1] != "base":
            continue
        ious = [float(x) for x in p[10:19]]
        rows[(p[4], p[5])] = {
            "z": float(p[2]),
            "area": float(p[6]),
            "b": float(p[7]),
            "iou": float(p[8]),
            "edge": int(float(p[9])),
            "ious": ious,
        }
    if len(rows) != 6:
        fail("expected 6 baseline curves in %s" % path.name)
    return rows


def fof_at(dx, dy, dz, mass, b_values):
    i_peak = sb.nucleus_index(dx, dy, dz)
    keeps = []
    for bpar in b_values:
        link = sb.link02(mass, float(bpar))
        root = sb.fof_largest(dx, dy, dz, link)
        keep = sb.peak_mask(root, dx, dy, dz, i_peak)
        keeps.append(keep)
        print("  b", float(bpar), "n", int(keep.sum()), flush=True)
    return keeps


def pack(triples):
    return " ".join(fmt(v) for triple in triples for v in triple)


def emit(fh, tag, z, nstar, los, band, triples):
    fh.write(
        "%s %.4f %d %s %s %s\n" % (tag, z, int(nstar), los, band, pack(triples))
    )


def curve(w, keeps, frac):
    return [scores_frac(w, keep, frac) for keep in keeps]


def orig_index():
    out = []
    for bpar in ORIG_B:
        k = int(np.argmin(np.abs(FINE_B - bpar)))
        if abs(float(FINE_B[k]) - float(bpar)) > 1e-8:
            fail("published node %s missing from the fine grid" % bpar)
        out.append(k)
    return out


def measure(step, nstar):
    cur = load_base(step, nstar)
    H = rv.prepare_halo(step, nstar)
    if int(H["nstar"]) != int(nstar):
        fail("catalogue nstar %s is not %s" % (H["nstar"], nstar))
    z = float(H["z"])
    nstar = int(H["nstar"])
    scale, n = rv.angular_scale(H["a_snap"], z)
    print("  arcsec scale", scale, "n", n, flush=True)
    kernels = []
    for amplitude in sb.WING_A:
        kern, _, _ = sb.psf_kernel(
            sb.PIX_AS, sb.FWHM_AS, sb.MOFFAT_BETA, amplitude, sb.KERNEL_AS)
        kernels.append(kern)
    raws, pixs = rv.project_raw(H["dx"], H["dy"], H["dz"], H["flux"], scale, n)
    images = []
    for ia in range(3):
        images.append([rv.convolve_band(raws[ia][ib], kernels[ib]) for ib in range(2)])
    del raws

    packed = []
    for ia, los in enumerate(rv.LOS):
        ix, iy, on = pixs[ia]
        for ib, band in enumerate(rv.BAND):
            row = cur[(los, band)]
            _mu, mask, labels, _nlab = rv.mu_labels(images[ia][ib], sb.THRESHOLDS[ib])
            cid, _fallback, _off = rv.component_at(labels, mask, n // 2, n // 2)
            region, area = rv.region_of(labels, cid, sb.PIX_AS)
            if abs(area - row["area"]) > 5e-4:
                fail("area %s %s got %.6f cur %.6f" % (los, band, area, row["area"]))
            inside = rv.sample_inside(region, ix, iy, on)
            if int(region.sum()) == 0:
                frac = np.zeros(H["flux"].shape[0], np.float64)
            else:
                frac_map = fftconvolve(region.astype(np.float64), kernels[ib], mode="same")
                frac = np.zeros(H["flux"].shape[0], np.float64)
                frac[on] = frac_map[iy[on], ix[on]]
                del frac_map
                frac[frac < 0.0] = 0.0
                if float(np.nanmax(frac)) > 1.0 + 1e-3:
                    fail("PSF fraction exceeds one on %s %s" % (los, band))
                w = H["flux"][:, ib]
                deposited = np.where(np.isfinite(w) & (w > 0.0), w, 0.0)
                pred = float(np.sum(deposited * frac))
                photo = float(images[ia][ib][region].sum())
                scale_p = max(abs(photo), abs(pred), 1e-30)
                if abs(photo - pred) > 1e-4 * scale_p:
                    fail("aperture sum %s %s photo %.8e pred %.8e" % (los, band, photo, pred))
            packed.append((los, band, row, inside, frac, ib))
            print("  iso", los, band, "area", area, flush=True)
    del images

    keeps = fof_at(H["dx"], H["dy"], H["dz"], H["mass"], FINE_B)
    oidx = orig_index()
    lines = []
    for los, band, row, inside, frac, ib in packed:
        w = H["flux"][:, ib]
        pos = [rv.light_scores(w, keeps[k], inside) for k in range(len(FINE_B))]
        psf = curve(w, keeps, frac)
        for j, k in enumerate(oidx):
            got = pos[k][0]
            ref = row["ious"][j]
            if np.isfinite(ref) != np.isfinite(got) or (np.isfinite(ref) and abs(got - ref) > 1e-4):
                fail("positional iou %s %s b %s got %s cur %s" % (los, band, ORIG_B[j], got, ref))
        lines.append(("POS", los, band, pos))
        lines.append(("PSF", los, band, psf))
        print("  scored", los, band, flush=True)
    del keeps

    rng = np.random.default_rng(np.random.SeedSequence([20261007, nstar]))
    draw = rng.random(H["dx"].size)
    thin_lines = []
    for tag, frac_keep in THIN:
        m = draw < frac_keep
        print("  thin", tag, "n", int(m.sum()), flush=True)
        tkeeps = fof_at(H["dx"][m], H["dy"][m], H["dz"][m], H["mass"][m], ORIG_B)
        for los, band, _row, inside, _frac, ib in packed:
            w = H["flux"][m, ib]
            inn = inside[m]
            triples = [rv.light_scores(w, keep, inn) for keep in tkeeps]
            thin_lines.append((tag, los, band, triples))
        del tkeeps

    outp = REV / f"fine_{step}_{nstar}.txt"
    tmp = REV / f"fine_{step}_{nstar}.txt.tmp"
    with tmp.open("w") as fh:
        fh.write("# Fine grid inserts 0.06 0.07 0.09 0.11 0.13 0.14. Ends remain 0.02 and 0.20.\n")
        fh.write("# POS is the positional score. PSF weights each star by the fraction of its kernel inside the adopted isophote.\n")
        fh.write("# THIN uses the published subset against that fixed isophote. The image is not rebuilt.\n")
        fh.write("# GRID_FINE %s\n" % " ".join("%.2f" % b for b in FINE_B))
        fh.write("# GRID_ORIG %s\n" % " ".join("%.2f" % b for b in ORIG_B))
        fh.write("# After the sightline: iou recall purity at each grid node, in grid order.\n")
        for tag, los, band, triples in lines:
            emit(fh, tag, z, nstar, los, band, triples)
        for tag, los, band, triples in thin_lines:
            emit(fh, tag, z, nstar, los, band, triples)
        fh.write("DONE fine %s %s\n" % (step, nstar))
    tmp.replace(outp)
    print("FINE_OK", step, nstar, flush=True)


def main():
    self_test()
    if len(sys.argv) == 2 and sys.argv[1] == "--self-test":
        return
    step = nstar = None
    args = sys.argv[1:]
    i = 0
    while i < len(args):
        if args[i] == "--step":
            step = args[i + 1]
            i += 2
        elif args[i] == "--nstar":
            nstar = args[i + 1]
            i += 2
        else:
            fail("unknown argument %s" % args[i])
    if step is None or nstar is None:
        fail("need --step and --nstar")
    outp = REV / f"fine_{step}_{nstar}.txt"
    if outp.is_file() and any(line.startswith("DONE ") for line in outp.read_text().splitlines()):
        print("skip", step, nstar, flush=True)
        return
    measure(step, nstar)


if __name__ == "__main__":
    main()
