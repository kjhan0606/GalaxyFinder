#!/usr/bin/env python3
"""Cross-validate the fixed linking parameter while holding out halo groups.

The exact-ID analysis covers outputs 00178--00786. A second, full-sample
sensitivity analysis adds reciprocal-nearest positional candidates at the
later outputs (and for the single z=3 halo). Those later assignments remain
provisional and are reported separately.
"""
from __future__ import annotations

import csv
import json
from collections import defaultdict, Counter
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
SUPP = ROOT / "supplement"
B_GRID = np.array([0.02, 0.03, 0.05, 0.08, 0.10, 0.12, 0.15, 0.18, 0.20])
BOX = 44.171267
NEIGHBOUR_LIMIT = 0.5
EXACT_INTERVALS = {("00178", "00475"), ("00475", "00786")}
POSITIONAL_INTERVALS = {("00786", "00967"), ("00967", "01001")}


def read_tsv(path: Path):
    with path.open() as stream:
        return list(csv.DictReader(
            (line for line in stream if line.strip() and not line.startswith("#")),
            delimiter="\t",
        ))


def write_tsv(path: Path, rows):
    if not rows:
        raise ValueError(f"No rows to write: {path}")
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter="\t",
                                lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def components(keys, edges):
    parent = {key: key for key in keys}

    def find(key):
        while parent[key] != key:
            parent[key] = parent[parent[key]]
            key = parent[key]
        return key

    for a, b, _ in edges:
        if a not in parent or b not in parent:
            continue
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[max(ra, rb)] = min(ra, rb)

    groups = defaultdict(list)
    for key in keys:
        groups[find(key)].append(key)
    component_of = {}
    for members in groups.values():
        label = min(members)
        label_string = f"{label[0]}_{label[1]}"
        for key in members:
            component_of[key] = label_string
    return component_of, groups


def periodic_distance(a, b):
    delta = np.asarray(a, dtype=float) - np.asarray(b, dtype=float)
    delta -= BOX * np.rint(delta / BOX)
    return float(np.linalg.norm(delta))


def make_edges(all_keys, positions):
    key_set = set(all_keys)
    exact_rows = read_tsv(SUPP / "lineage_dm_id_overlaps.tsv")
    exact_edges = []
    for row in exact_rows:
        interval = (row["earlier_step"], row["later_step"])
        if interval not in EXACT_INTERVALS:
            continue
        if int(row["is_earlier_best_descendant"]) != 1:
            continue
        if float(row["fraction_earlier_dm_in_later"]) <= 0.5:
            continue
        a = (row["earlier_step"], int(row["earlier_nstar"]))
        b = (row["later_step"], int(row["later_nstar"]))
        if a in key_set and b in key_set:
            exact_edges.append((a, b, "exact_DM_particle_ID"))

    positional_edges = []
    for row in read_tsv(SUPP / "lineage_position_candidates.tsv"):
        interval = (row["earlier_step"], row["later_step"])
        if interval not in POSITIONAL_INTERVALS or int(row["provisional_link"]) != 1:
            continue
        a = (row["earlier_step"], int(row["earlier_nstar"]))
        b = (row["later_step"], int(row["later_nstar"]))
        if a in key_set and b in key_set:
            positional_edges.append((a, b, "provisional_reciprocal_position"))

    # The artificial-redshift table starts at z=1.99. Add its omitted z=3
    # entry and use the same reciprocal-nearest, <3 cMpc/h positional rule.
    z3 = [key for key in all_keys if key[0] == "00081"]
    z199 = [key for key in all_keys if key[0] == "00178"]
    if z3 and z199:
        distances = np.asarray([[periodic_distance(positions[a], positions[b])
                                 for b in z199] for a in z3])
        reverse_nearest = np.argmin(distances, axis=0)
        for i, a in enumerate(z3):
            j = int(np.argmin(distances[i]))
            if distances[i, j] < 3.0 and int(reverse_nearest[j]) == i:
                positional_edges.append((a, z199[j], "provisional_reciprocal_position"))

    return exact_edges, positional_edges


def make_curves(keys):
    rows_by_view = defaultdict(list)
    for row in read_tsv(SUPP / "recomputed_curves.tsv"):
        if row["variant"] != "consistent_distance":
            continue
        key = (row["step"], int(row["nstar"]), row["los"], row["band"])
        rows_by_view[key].append(row)

    curves = {}
    for step, nstar in keys:
        for band in ("V", "g"):
            view_curves = []
            for los in "xyz":
                rows = sorted(rows_by_view[(step, nstar, los, band)],
                              key=lambda row: float(row["b"]))
                if len(rows) != len(B_GRID):
                    raise ValueError(f"Incomplete curve for {step}:{nstar} {los} {band}")
                if not np.allclose([float(row["b"]) for row in rows], B_GRID):
                    raise ValueError(f"Unexpected b grid for {step}:{nstar} {los} {band}")
                view_curves.append(np.asarray([float(row["iou"]) for row in rows]))
            curves[(step, nstar, band)] = np.median(view_curves, axis=0)
    return curves


def run_group_cv(name, data_keys, component_of, group_members, retained, positions,
                 curves, evidence_label):
    held_keys = [key for key in data_keys if retained[key]]
    components_in_sample = sorted({component_of[key] for key in held_keys})
    keys_by_component = defaultdict(list)
    for key in data_keys:
        keys_by_component[component_of[key]].append(key)
    trainable_by_component = defaultdict(list)
    for key in held_keys:
        trainable_by_component[component_of[key]].append(key)

    folds, predictions, component_rows = [], [], []
    for key in data_keys:
        component_rows.append({
            "validation_definition": name,
            "step": key[0], "nstar": key[1],
            "z": positions["z"][key],
            "calibration_retained": int(retained[key]),
            "component": component_of[key],
            "link_evidence": evidence_label.get(component_of[key], "no_link"),
        })

    for held_component in components_in_sample:
        held = trainable_by_component[held_component]
        held_component_keys = keys_by_component[held_component]
        held_positions = [positions["xyz"][key] for key in held_component_keys]
        train_components, neighbour_removed_components = [], []
        for candidate in components_in_sample:
            if candidate == held_component:
                continue
            candidate_positions = [positions["xyz"][key]
                                   for key in keys_by_component[candidate]]
            close = min(periodic_distance(x, y)
                        for x in held_positions for y in candidate_positions) < NEIGHBOUR_LIMIT
            if close:
                neighbour_removed_components.append(candidate)
            else:
                train_components.append(candidate)
        training = [key for comp in train_components
                    for key in trainable_by_component[comp]]
        if not train_components or not training:
            raise ValueError(f"No training data in {name} fold {held_component}")

        for band in ("V", "g"):
            equal_component_curves = [
                np.median([curves[(*key, band)] for key in trainable_by_component[comp]], axis=0)
                for comp in train_components
            ]
            objective = np.median(equal_component_curves, axis=0)
            selected_index = int(np.argmax(objective))
            selected_b = float(B_GRID[selected_index])
            folds.append({
                "validation_definition": name,
                "held_out_component": held_component,
                "band": band,
                "n_train_entries": len(training),
                "n_train_components": len(train_components),
                "n_neighbour_components_removed": len(neighbour_removed_components),
                "n_neighbour_entries_removed": sum(len(trainable_by_component[c])
                                                    for c in neighbour_removed_components),
                "n_held_entries": len(held),
                "selected_b": selected_b,
            })
            for key in held:
                curve = curves[(*key, band)]
                predictions.append({
                    "validation_definition": name,
                    "component": held_component,
                    "step": key[0], "nstar": key[1],
                    "z": positions["z"][key],
                    "band": band,
                    "selected_b": selected_b,
                    "held_out_iou": float(curve[selected_index]),
                    "fixed_b010_iou": float(curve[int(np.flatnonzero(B_GRID == 0.10)[0])]),
                    "fixed_b020_iou": float(curve[-1]),
                })

    summaries = {}
    for band_index, band in enumerate(("V", "g")):
        rows = [row for row in predictions if row["band"] == band]
        by_component = defaultdict(list)
        for row in rows:
            by_component[row["component"]].append(row)
        group_stats = []
        for component, group_rows in by_component.items():
            group_stats.append({
                "component": component,
                "selected_iou": float(np.median([r["held_out_iou"] for r in group_rows])),
                "b010_iou": float(np.median([r["fixed_b010_iou"] for r in group_rows])),
                "b020_iou": float(np.median([r["fixed_b020_iou"] for r in group_rows])),
            })
        gains_selected = np.asarray([r["selected_iou"] - r["b020_iou"] for r in group_stats])
        gains_b010 = np.asarray([r["b010_iou"] - r["b020_iou"] for r in group_stats])
        rng = np.random.default_rng(20261030 + band_index + (0 if name == "exact_verified" else 100))

        def bootstrap_interval(values):
            draws = rng.integers(0, len(values), size=(5000, len(values)))
            means = values[draws].mean(axis=1)
            return [float(x) for x in np.quantile(means, [0.16, 0.84])]

        selected_values = [float(x) for x in sorted({r["selected_b"] for r in folds
                                                     if r["band"] == band})]
        counts = Counter(str(r["selected_b"]) for r in folds if r["band"] == band)
        summaries[band] = {
            "n_retained_entries": len(rows),
            "n_held_out_components": len(group_stats),
            "selected_b_values": selected_values,
            "selected_b_fold_counts": dict(counts),
            "median_held_out_iou_equal_component": float(np.median(
                [r["selected_iou"] for r in group_stats])),
            "mean_held_out_iou_equal_component": float(np.mean(
                [r["selected_iou"] for r in group_stats])),
            "mean_gain_selected_b_over_020_equal_component": float(np.mean(gains_selected)),
            "gain_selected_b_over_020_bootstrap_16_84": bootstrap_interval(gains_selected),
            "mean_gain_fixed_010_over_020_equal_component": float(np.mean(gains_b010)),
            "gain_fixed_010_over_020_bootstrap_16_84": bootstrap_interval(gains_b010),
        }
    return folds, predictions, component_rows, summaries


def main():
    manifest = read_tsv(SUPP / "sample_manifest.tsv")
    all_keys = [(row["step"], int(row["nstar"])) for row in manifest]
    retained = {(row["step"], int(row["nstar"])):
                row["calibration_retained"].strip().lower() in {"true", "1"}
                for row in manifest}
    position_rows = read_tsv(SUPP / "halo_positions.tsv")
    positions = {
        "xyz": {(r["step"], int(r["nstar"])):
                np.asarray([float(r[k]) for k in
                            ("x_hinv_cMpc", "y_hinv_cMpc", "z_hinv_cMpc")])
                for r in position_rows},
        "z": {(r["step"], int(r["nstar"])): float(r["z"])
              for r in manifest},
    }
    if set(all_keys) != set(positions["xyz"]) or len(all_keys) != 76:
        raise ValueError("Expected one position for each of the 76 manifest entries")
    exact_edges, positional_edges = make_edges(all_keys, positions["xyz"])
    curves = make_curves(all_keys)

    exact_keys = [key for key in all_keys if key[0] in {"00178", "00475", "00786"}]
    exact_component, exact_groups = components(exact_keys, exact_edges)
    exact_evidence = {}
    for a, _, evidence in exact_edges:
        exact_evidence[exact_component[a]] = "exact_DM_particle_ID"
    exact_retained = {key: retained[key] for key in exact_keys}
    exact_result = run_group_cv(
        "exact_verified", exact_keys, exact_component, exact_groups,
        exact_retained, positions, curves, exact_evidence,
    )

    hybrid_edges = exact_edges + positional_edges
    hybrid_component, hybrid_groups = components(all_keys, hybrid_edges)
    evidence_by_component = defaultdict(set)
    for a, _, evidence in hybrid_edges:
        evidence_by_component[hybrid_component[a]].add(evidence)
    hybrid_evidence = {}
    for component in set(hybrid_component.values()):
        evidence = evidence_by_component[component]
        if "exact_DM_particle_ID" in evidence and "provisional_reciprocal_position" in evidence:
            hybrid_evidence[component] = "exact+provisional_position"
        elif "exact_DM_particle_ID" in evidence:
            hybrid_evidence[component] = "exact_DM_particle_ID"
        elif "provisional_reciprocal_position" in evidence:
            hybrid_evidence[component] = "provisional_reciprocal_position"
        else:
            hybrid_evidence[component] = "no_link"
    hybrid_result = run_group_cv(
        "hybrid_provisional_full", all_keys, hybrid_component, hybrid_groups,
        retained, positions, curves, hybrid_evidence,
    )

    fold_rows = exact_result[0] + hybrid_result[0]
    prediction_rows = exact_result[1] + hybrid_result[1]
    component_rows = exact_result[2] + hybrid_result[2]
    summary = {
        "description": "The training objective is the median overlap curve with equal weight per component; each fold holds out an entire component and removes any training component with a member within 0.5 h^-1 cMpc of a held-out member.",
        "exact_verified_scope": {
            "outputs": ["00178", "00475", "00786"],
            "entries": len(exact_keys),
            "accepted_exact_DM_ID_edges": len(exact_edges),
            "components": len(exact_groups),
            "position_links_used": 0,
            "results_by_band": exact_result[3],
        },
        "hybrid_full_scope": {
            "entries": len(all_keys),
            "calibration_retained_entries": sum(retained.values()),
            "accepted_exact_DM_ID_edges": len(exact_edges),
            "accepted_provisional_position_edges": len(positional_edges),
            "components": len(hybrid_groups),
            "components_with_positional_candidates": sum(
                evidence.startswith("provisional") or evidence == "exact+provisional_position"
                for evidence in hybrid_evidence.values()),
            "results_by_band": hybrid_result[3],
        },
    }
    write_tsv(SUPP / "lineage_cv_folds.tsv", fold_rows)
    write_tsv(SUPP / "lineage_cv_predictions.tsv", prediction_rows)
    write_tsv(SUPP / "lineage_cv_components.tsv", component_rows)
    (SUPP / "lineage_cv_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
