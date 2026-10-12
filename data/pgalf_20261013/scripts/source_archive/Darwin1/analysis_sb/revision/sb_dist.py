#!/usr/bin/env python3
"""Residual light versus signed distance from the adopted isophote.

Does not rewrite rev_*.txt and does not refit the linking length.
The isophote is the baseline adopted component. The group is the
friends-of-friends set at that sightline's baseline grid optimum.
Positive distance is outside the isophote. The score is the whole-star
flux of a star whose pixel lies at that distance. An empty isophote has
no boundary.
"""
import sys
from pathlib import Path

import numpy as np
from scipy import ndimage

REV = Path("/gpfs/kjhan/Others/Jaehyun/2026.Apr/Darwin1/analysis_sb/revision")
sys.path.insert(0, str(REV))
import sb_revise as rv
import sb_vs_b as sb


def fail(msg):
    print("DIST_FAIL", msg, flush=True)
    raise SystemExit(2)


def signed_arcsec(region):
    if int(region.sum()) == 0:
        return None
    inside = ndimage.distance_transform_edt(region)
    outside = ndimage.distance_transform_edt(~region)
    return (outside - inside) * float(sb.PIX_AS)


def n_border(region):
    if region.size == 0 or not np.any(region):
        return 0
    edge = np.zeros(region.shape, dtype=bool)
    edge[0, :] = True
    edge[-1, :] = True
    edge[:, 0] = True
    edge[:, -1] = True
    return int(np.count_nonzero(region & edge))


def light(w, mask):
    if not np.any(mask):
        return 0.0
    return float(w[mask].sum())


def leq(a, b):
    scale = max(abs(a), abs(b), 1e-30)
    return a <= b + 1e-6 * scale


def self_test():
    region = np.zeros((7, 7), dtype=bool)
    region[2:5, 2:5] = True
    d = signed_arcsec(region)
    pix = float(sb.PIX_AS)
    if abs(d[2, 3] + pix) > 1e-9:
        fail("interior edge distance %s" % d[2, 3])
    if abs(d[1, 3] - pix) > 1e-9:
        fail("exterior face distance %s" % d[1, 3])
    if abs(d[1, 1] - pix * np.sqrt(2.0)) > 1e-9:
        fail("exterior diagonal distance %s" % d[1, 1])
    if abs(d[3, 3] + 2.0 * pix) > 1e-9:
        fail("interior centre distance %s" % d[3, 3])
    if signed_arcsec(np.zeros((5, 5), dtype=bool)) is not None:
        fail("empty isophote returned a distance")
    ring = rv.ring_of(region)
    ext = ring & ~region
    inn = ring & region
    if not np.all((d[ext] > 0.0) & (d[ext] <= 1.0)):
        fail("exterior ring is outside 1 arcsec")
    if not np.all((d[inn] < 0.0) & (d[inn] >= -1.0)):
        fail("interior ring is outside 1 arcsec")
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
        if len(p) < 10 or p[0] != "CUR" or p[1] != "base":
            continue
        los, band = p[4], p[5]
        rows[(los, band)] = {
            "z": float(p[2]),
            "area": float(p[6]),
            "b": float(p[7]),
            "iou": float(p[8]),
            "edge": int(float(p[9])),
        }
    if len(rows) != 6:
        fail("expected 6 baseline curves in %s" % path.name)
    return rows


def snap_index(b):
    k = int(np.argmin(np.abs(rv.B_GRID - b)))
    if abs(float(rv.B_GRID[k]) - float(b)) > 1e-4:
        fail("b %s is off the grid" % b)
    return k


def fof_selected(dx, dy, dz, mass, indices):
    i_peak = sb.nucleus_index(dx, dy, dz)
    keeps = {}
    for k in indices:
        bpar = float(rv.B_GRID[k])
        link = sb.link02(mass, bpar)
        root = sb.fof_largest(dx, dy, dz, link)
        keeps[k] = sb.peak_mask(root, dx, dy, dz, i_peak)
        print("  b", bpar, "n", int(keeps[k].sum()), flush=True)
    return keeps


def fmt(x):
    if not np.isfinite(x):
        return "nan"
    return "%.8e" % float(x)


def emit(fh, cur, los, band, area, ntouch, empty, got, cols):
    fh.write(
        "DIST %s %s %s %s %.3f %.4f %d %s %s %d %d %s\n"
        % (
            "%.4f" % cur["z"],
            cur_nstar(cur),
            los,
            band,
            area,
            float(rv.B_GRID[snap_index(cur["b"])]),
            cur["edge"],
            fmt(cur["iou"]),
            fmt(got),
            ntouch,
            empty,
            " ".join(fmt(c) for c in cols),
        )
    )


def cur_nstar(cur):
    return cur["nstar"]


def measure(step, nstar):
    cur = load_base(step, nstar)
    for row in cur.values():
        row["nstar"] = int(nstar)
    H = rv.prepare_halo(step, nstar)
    z = H["z"]
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
    indices = set()
    for ia, los in enumerate(rv.LOS):
        ix, iy, on = pixs[ia]
        for ib, band in enumerate(rv.BAND):
            row = cur[(los, band)]
            mu, mask, labels, _nlab = rv.mu_labels(images[ia][ib], sb.THRESHOLDS[ib])
            cid, _fallback, _off = rv.component_at(labels, mask, n // 2, n // 2)
            region, area = rv.region_of(labels, cid, sb.PIX_AS)
            if abs(area - row["area"]) > 5e-4:
                fail(
                    "area %s %s got %.6f cur %.6f" % (los, band, area, row["area"])
                )
            k = snap_index(row["b"])
            indices.add(k)
            packed.append((ia, ib, los, band, row, region, area, k, ix, iy, on))
            print("  iso", los, band, "area", area, "border", n_border(region), flush=True)
    del images

    keeps = fof_selected(H["dx"], H["dy"], H["dz"], H["mass"], sorted(indices))
    lines = []
    for ia, ib, los, band, row, region, area, k, ix, iy, on in packed:
        w = H["flux"][:, ib]
        keep = keeps[k]
        inside = rv.sample_inside(region, ix, iy, on)
        touch = n_border(region)
        empty = int(area <= 0.0 or not np.any(region))
        if empty:
            grp = light(w, keep)
            cols = [grp, 0.0, 0.0, grp, 0.0] + [np.nan] * 9
            lines.append((row, los, band, area, touch, 1, np.nan, cols))
            print("  empty", los, band, flush=True)
            continue
        signed = signed_arcsec(region)
        d = np.full(w.shape[0], np.nan)
        d[on] = signed[iy[on], ix[on]]
        del signed
        if np.any(inside & on & ~(d < 0.0)):
            fail("inside star has a non-negative distance %s %s" % (los, band))
        if np.any((~inside) & on & ~(d > 0.0)):
            fail("outside star has a non-positive distance %s %s" % (los, band))
        iou, _rec, _pur = rv.light_scores(w, keep, inside)
        if not np.isfinite(iou) or abs(iou - row["iou"]) > 1e-4:
            fail("iou %s %s got %s cur %s" % (los, band, iou, row["iou"]))
        ring = rv.sample_inside(rv.ring_of(region), ix, iy, on)
        miss = keep & ~inside
        extra = inside & ~keep
        on_miss = miss & on
        grp = light(w, keep)
        iso = light(w, inside)
        both = light(w, keep & inside)
        lmiss = light(w, miss)
        lextra = light(w, extra)
        miss_ring = light(w, miss & ring)
        miss_1 = light(w, on_miss & (d <= 1.0))
        miss_2 = light(w, on_miss & (d <= 2.0))
        miss_5 = light(w, on_miss & (d <= 5.0))
        miss_off = light(w, miss & ~on)
        extra_ring = light(w, extra & ring)
        extra_1 = light(w, extra & (d >= -1.0))
        extra_2 = light(w, extra & (d >= -2.0))
        extra_5 = light(w, extra & (d >= -5.0))
        cols = [
            grp, iso, both, lmiss, lextra,
            miss_ring, miss_1, miss_2, miss_5, miss_off,
            extra_ring, extra_1, extra_2, extra_5,
        ]
        if not leq(miss_ring, miss_1) or not leq(miss_1, miss_2) or not leq(miss_2, miss_5):
            fail("miss distances are not nested %s %s" % (los, band))
        if not leq(miss_5 + miss_off, lmiss):
            fail("miss split exceeds the outside group light %s %s" % (los, band))
        if not leq(extra_ring, extra_1) or not leq(extra_1, extra_2) or not leq(extra_2, extra_5):
            fail("extra distances are not nested %s %s" % (los, band))
        if not leq(extra_5, lextra):
            fail("extra split exceeds the extra isophote light %s %s" % (los, band))
        lines.append((row, los, band, area, touch, 0, iou, cols))
        print("  dist", los, band, "iou", iou, "miss", lmiss, "extra", lextra, flush=True)

    outp = REV / f"dist_{step}_{nstar}.txt"
    tmp = REV / f"dist_{step}_{nstar}.txt.tmp"
    with tmp.open("w") as fh:
        fh.write("# Signed distance from the adopted isophote. Positive is outside.\n")
        fh.write("# b is the baseline optimum. It is not refitted. edge 1 is not a matched length.\n")
        fh.write("# empty 1 has no boundary, and the distance columns are nan.\n")
        fh.write("# Columns after empty: L_grp L_iso L_both L_miss L_extra\n")
        fh.write("# L_miss_ring L_miss_le1 L_miss_le2 L_miss_le5 L_miss_off\n")
        fh.write("# L_extra_ring L_extra_le1 L_extra_le2 L_extra_le5\n")
        fh.write("# le1, le2, and le5 are Euclidean arcsec on that side of the boundary.\n")
        fh.write("# ring is the 8-connected morphological ring. L_miss_off is off the stamp.\n")
        for row, los, band, area, touch, empty, got, cols in lines:
            emit(fh, row, los, band, area, touch, empty, got, cols)
        fh.write("DONE dist %s %s\n" % (step, int(nstar)))
    tmp.replace(outp)
    print("DIST_OK", step, nstar, flush=True)


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
    outp = REV / f"dist_{step}_{nstar}.txt"
    if outp.is_file() and any(line.startswith("DONE ") for line in outp.read_text().splitlines()):
        print("skip", step, nstar, flush=True)
        return
    measure(step, nstar)


if __name__ == "__main__":
    main()
