#!/usr/bin/env python3
"""Summarize the corrected boundary and distance calculations from saved sums.

The original sample flags and per-view selected parameters remain authoritative
for baseline comparisons. New distance-corrected parameters use unrounded scores.
No particle files or friends-of-friends library are needed.
"""
import csv
import json
from collections import defaultdict
from pathlib import Path

import numpy as np

SUPP = Path(__file__).resolve().parents[1] / "supplement"


def read(name):
    return list(csv.DictReader(
        [line for line in (SUPP / name).read_text().splitlines()
         if line.strip() and not line.startswith("#")], delimiter="\t"))


def quantiles(values):
    values = np.asarray(values, dtype=float)
    finite = values[np.isfinite(values)]
    return dict(n=len(values), n_finite=len(finite),
                p16_median_p84=np.quantile(finite, [.16, .5, .84]).tolist()
                if len(finite) else [])


def write(name, rows):
    with (SUPP / name).open("w") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main():
    metadata = read("halo_metadata.tsv")
    good = {(r["step"], int(r["nstar"])):
            bool(int(r["interior_V"]) and int(r["interior_g"])) for r in metadata}
    originals = {(r["step"], int(r["nstar"]), r["los"], r["band"]): r
                 for r in read("boundary_ring.tsv")}
    curves = defaultdict(list)
    for row in read("recomputed_curves.tsv"):
        key = row["step"], int(row["nstar"]), row["los"], row["band"], row["variant"]
        curves[key].append(row)
    if len(originals) != 6 * len(metadata) or len(curves) != 12 * len(metadata):
        raise ValueError("The revision summaries require all 76 halo entries.")

    comparisons = []
    for key, original in originals.items():
        old_b = float(original["b_best"])
        baseline = next(r for r in curves[(*key, "baseline")]
                        if np.isclose(float(r["b"]), old_b))
        changed = sorted(curves[(*key, "consistent_distance")],
                         key=lambda r: float(r["b"]))
        scores = np.array([float(r["iou"]) for r in changed])
        finite = np.isfinite(scores)
        imax = int(np.argmax(np.where(finite, scores, -np.inf)))
        best = changed[imax]
        area = float(baseline["area_arcsec2"])
        new_area = float(best["area_arcsec2"])
        comparisons.append(dict(
            step=key[0], nstar=key[1], los=key[2], band=key[3],
            baseline_halo_interior=int(good[key[:2]]),
            baseline_view_edge=int(original["edge"]),
            corrected_view_edge=int(not finite.any() or imax in (0, len(changed)-1)),
            baseline_b=old_b, corrected_b=float(best["b"]),
            delta_b=float(best["b"])-old_b,
            baseline_iou=float(baseline["iou"]), corrected_iou=float(best["iou"]),
            delta_iou=float(best["iou"])-float(baseline["iou"]),
            baseline_area_arcsec2=area, corrected_area_arcsec2=new_area,
            area_ratio=new_area/area if area else float("nan")))
    write("distance_sensitivity.tsv", comparisons)

    summary = {}
    for band in ("V", "g"):
        halos = defaultdict(list)
        for row in comparisons:
            if row["band"] == band:
                halos[row["step"], row["nstar"]].append(row)
        paired = []
        retained = []
        lost = []
        gained = []
        for key, rows in halos.items():
            current = len(rows) == 3 and not any(r["corrected_view_edge"] for r in rows)
            if current:
                retained.append(np.median([r["corrected_b"] for r in rows]))
            if good[key] and not current:
                lost.append(list(key))
            if current and not good[key]:
                gained.append(list(key))
            if not (current and good[key]):
                continue
            old_b = np.median([r["baseline_b"] for r in rows])
            new_b = np.median([r["corrected_b"] for r in rows])
            old_iou = np.median([r["baseline_iou"] for r in rows])
            new_iou = np.median([r["corrected_iou"] for r in rows])
            paired.append((old_b, new_b, old_iou, new_iou))
        paired = np.asarray(paired)
        halo_rings = defaultdict(list)
        for key, row in originals.items():
            if key[3] == band and good[key[:2]]:
                halo_rings[key[:2]].append(row)
        boundary = {}
        for node in ("best", "0.10", "0.20"):
            for measure in ("iou", "completeness", "precision"):
                column = f"ring_{measure}_{node}"
                boundary[column] = quantiles([
                    np.median([float(row[column]) for row in rows])
                    for rows in halo_rings.values()])
        views = [row for row in comparisons if row["band"] == band]
        per_output = {}
        for step in sorted({key[0] for key in halos}):
            vals = [np.median([r["corrected_b"] for r in rows])
                    for key, rows in halos.items() if key[0] == step
                    and not any(r["corrected_view_edge"] for r in rows)]
            per_output[step] = quantiles(vals)
        summary[band] = dict(
            boundary=boundary,
            distance=dict(
                corrected_halo_b=quantiles(retained), per_output_b=per_output,
                n_paired=len(paired),
                n_paired_b_changed=int(np.count_nonzero(abs(paired[:, 1]-paired[:, 0]) > 1e-8)),
                paired_delta_b=quantiles(paired[:, 1]-paired[:, 0]),
                paired_delta_iou=quantiles(paired[:, 3]-paired[:, 2]),
                corrected_paired_iou=quantiles(paired[:, 3]),
                originally_retained_lost=lost, newly_retained=gained,
                n_empty_views=sum(r["corrected_area_arcsec2"] == 0 for r in views),
                n_area_increases=sum(r["corrected_area_arcsec2"] > r["baseline_area_arcsec2"] + 1e-8
                                     for r in views),
                view_area_ratio=quantiles([r["area_ratio"] for r in views])))

    catalogue = read("catalogue_audit.tsv")
    matches = lambda name: sum(int(r["selected_index"]) == int(r[name]) for r in catalogue)
    summary["catalogue"] = dict(
        n_halos=len(catalogue),
        reconstructed_count_matches=sum(int(r["nstar"]) == int(r["reconstructed_count"]) for r in catalogue),
        maximum_relative_mass_difference=max(abs(float(r["summed_stellar_mass"])
                                                 / float(r["catalogue_stellar_mass"])-1)
                                             for r in catalogue),
        central_is_brightest_V=matches("brightest_V"),
        central_is_brightest_g=matches("brightest_g"),
        central_is_most_massive=matches("most_massive"))
    (SUPP / "revision_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
