#!/usr/bin/env python3
"""Fetch exact FoF stellar membership for one SubLink progenitor epoch.

For each host, this saves the parent FoF star cutout and a boolean mask that
removes every non-central subhalo with at least 100 stellar particles. The
central is fixed by the SubLink branch's GroupFirstSub, not re-selected from
the particle data. Particle IDs are checked exactly against the parent cutout.

The TNG API key is read from ~/.config/illustris-tng/api_key and is never
written to output or included in diagnostics.
"""
from __future__ import annotations

import argparse
import csv
import io
import json
import os
import tempfile
import time
from pathlib import Path

import h5py
import numpy as np
import requests

HERE = Path(__file__).resolve().parent
LINEAGE = HERE / "lineage_epochs.csv"
COHORT = HERE / "lineage_cohort.csv"
KEYFILE = Path.home() / ".config/illustris-tng/api_key"
OUT = HERE / "tree_epoch_membership"
STAR_FIELDS = (
    "Coordinates,Masses,ParticleIDs,GFM_InitialMass,GFM_Metallicity,"
    "GFM_StellarFormationTime"
)


def safe_get(session: requests.Session, url: str, params=None, timeout=(20, 180)):
    last_status = None
    for attempt in range(6):
        try:
            response = session.get(url, params=params, timeout=timeout)
            if response.status_code in (403, 429, 500, 502, 503, 504):
                last_status = f"HTTP {response.status_code}"
                time.sleep(min(2 ** attempt, 30))
                continue
            if not response.ok:
                raise RuntimeError(f"TNG API returned HTTP {response.status_code}")
            return response
        except requests.RequestException as exc:
            # Avoid exception strings: redirected URLs may carry temporary credentials.
            last_status = type(exc).__name__
            if attempt == 5:
                break
            time.sleep(min(2 ** attempt, 30))
    raise RuntimeError(f"TNG request failed after retries ({last_status})")


def read_rows(path: Path):
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream))


def epoch_row(snapshot: int, rank: int):
    rows = [r for r in read_rows(LINEAGE)
            if int(r["snapshot"]) == snapshot and int(r["sample_rank"]) == rank]
    if len(rows) != 1 or int(rows[0].get("branch_node_present", 0)) != 1:
        raise RuntimeError(f"No unique SubLink node for rank={rank}, snapshot={snapshot}")
    row = rows[0]
    if int(row["is_fof_central"]) != 1:
        raise RuntimeError(f"Rank={rank} is not an FoF central at snapshot={snapshot}")
    return row


def fetch_parent(session, api_root: str, snapshot: int, group_id: int,
                 snapshot84_parent_dir: Path | None = None):
    if snapshot == 84:
        if snapshot84_parent_dir is None:
            raise RuntimeError(
                "Snapshot 84 retrieval needs the original parent membership cutout. "
                "Supply --snapshot84-parent-dir; the derived snapshot-84 membership "
                "used in this paper is already included in tree_epoch_membership/.")
        source = snapshot84_parent_dir / f"halo_{group_id:08d}_stars_membership.h5"
        if not source.is_file():
            raise FileNotFoundError(source)
        with h5py.File(source, "r") as handle:
            if int(handle.attrs["GroupID"]) != group_id or int(handle.attrs["Snapshot"]) != 84:
                raise RuntimeError(f"Cached descendant cutout metadata mismatch: {source.name}")
            stars = {key: handle["PartType4"][key][:] for key in (
                "Coordinates", "Masses", "ParticleIDs", "GFM_InitialMass",
                "GFM_Metallicity", "GFM_StellarFormationTime")}
            old_remove_satellite = handle["PartType4/RemoveSatelliteGalaxy"][:].astype(bool)
            old_central_id = int(handle.attrs["CentralSubhaloID"])
            old_central_nstar = int(handle.attrs["CentralNstar"])
            old_parent_bytes = int(handle.attrs["ParentCutoutBytes"])
            header = {"Redshift": float(handle.attrs["Redshift"]),
                      "BoxSize": float(handle.attrs.get("BoxSize", 205000.0)),
                      "CutoutID": group_id}
        return stars, header, old_parent_bytes, old_remove_satellite, old_central_id, old_central_nstar
    response = safe_get(session, f"{api_root}halos/{group_id}/cutout.hdf5",
                        params={"stars": STAR_FIELDS})
    with h5py.File(io.BytesIO(response.content), "r") as handle:
        header = {key: handle["Header"].attrs[key] for key in handle["Header"].attrs}
        if "PartType4" not in handle:
            raise RuntimeError(f"FoF group {group_id}: no PartType4 in cutout")
        stars = {key: handle["PartType4"][key][:] for key in (
            "Coordinates", "Masses", "ParticleIDs", "GFM_InitialMass",
            "GFM_Metallicity", "GFM_StellarFormationTime")}
    ids = stars["ParticleIDs"]
    if not len(ids) or len(np.unique(ids)) != len(ids):
        raise RuntimeError(f"FoF group {group_id}: empty or duplicate star ParticleIDs")
    if any(len(value) != len(ids) for value in stars.values()):
        raise RuntimeError(f"FoF group {group_id}: inconsistent star-array lengths")
    cutout_id = int(header.get("CutoutID", group_id))
    if cutout_id != group_id:
        raise RuntimeError(f"FoF cutout ID mismatch: requested {group_id}, got {cutout_id}")
    return stars, header, len(response.content), None, None, None


def candidate_subhalos(session, api_root: str, group_id: int):
    response = safe_get(session, f"{api_root}subhalos/", params={
        "grnr": group_id, "len_stars__gte": 100, "limit": 1000
    })
    page = response.json()
    rows = page.get("results", [])
    if page.get("next") is not None or page.get("count") != len(rows):
        raise RuntimeError(
            f"FoF group {group_id}: incomplete >=100-star subhalo list "
            f"({page.get('count')} listed, {len(rows)} returned)"
        )
    ids = [int(item["id"]) for item in rows]
    if not ids or len(ids) != len(set(ids)):
        raise RuntimeError(f"FoF group {group_id}: empty or duplicate candidate list")
    return ids


def subhalo_star_ids(session, api_root: str, subhalo_id: int):
    response = safe_get(session, f"{api_root}subhalos/{subhalo_id}/cutout.hdf5",
                        params={"stars": "ParticleIDs"})
    with h5py.File(io.BytesIO(response.content), "r") as handle:
        if "PartType4/ParticleIDs" not in handle:
            raise RuntimeError(f"Subhalo {subhalo_id}: no stellar ParticleIDs")
        ids = handle["PartType4/ParticleIDs"][:]
        cutout_id = int(handle["Header"].attrs.get("CutoutID", -1))
    if cutout_id != subhalo_id or len(ids) < 100 or len(np.unique(ids)) != len(ids):
        raise RuntimeError(f"Subhalo {subhalo_id}: invalid ID-only cutout")
    return ids


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--snapshot", type=int, choices=(84, 67, 50), required=True)
    parser.add_argument("--sample-rank", type=int, required=True)
    parser.add_argument("--snapshot84-parent-dir", type=Path, default=None,
                        help="Directory containing original snapshot-84 FoF membership cutouts")
    args = parser.parse_args()
    if not KEYFILE.is_file() or not KEYFILE.read_text().strip():
        raise RuntimeError("Protected TNG API key is missing or empty")
    if not COHORT.is_file() or not LINEAGE.is_file():
        raise FileNotFoundError("Extract SubLink epoch tables before fetching cutouts")
    selected = {int(r["sample_rank"]) for r in read_rows(COHORT)
                if int(r["central_at_all_targets"]) == 1}
    if args.sample_rank not in selected:
        raise RuntimeError(f"Sample rank {args.sample_rank} is not in the locked 36-tree cohort")

    epoch = epoch_row(args.snapshot, args.sample_rank)
    group_id = int(epoch["progenitor_fof_group_id"])
    central_id = int(epoch["progenitor_subhalo_id"])
    outpath = OUT / f"snap{args.snapshot:03d}_rank{args.sample_rank:02d}_group{group_id:06d}.h5"
    if outpath.exists():
        raise FileExistsError(f"Refusing to overwrite existing output: {outpath}")
    OUT.mkdir(parents=True, exist_ok=True)

    key = KEYFILE.read_text().strip()
    session = requests.Session()
    session.headers.update({"api-key": key, "Accept": "application/json,application/octet-stream"})
    del key
    api_root = f"https://www.tng-project.org/api/TNG300-1/snapshots/{args.snapshot}/"

    stars, header, parent_bytes, old_remove_satellite, old_central_id, old_central_nstar = \
        fetch_parent(session, api_root, args.snapshot, group_id,
                     args.snapshot84_parent_dir.resolve() if args.snapshot84_parent_dir else None)
    parent_ids = stars["ParticleIDs"]
    order = np.argsort(parent_ids)
    sorted_parent_ids = parent_ids[order]
    details = []
    central_parent = None
    if args.snapshot == 84:
        if old_central_id != central_id or old_remove_satellite.size != parent_ids.size:
            raise RuntimeError(f"Cached descendant membership does not match SubLink central {central_id}")
        ids = subhalo_star_ids(session, api_root, central_id)
        locations = np.searchsorted(sorted_parent_ids, ids)
        if np.any(locations >= sorted_parent_ids.size):
            raise RuntimeError(f"Subhalo {central_id}: IDs are not all in parent FoF cutout")
        if not np.array_equal(sorted_parent_ids[locations], ids):
            raise RuntimeError(f"Subhalo {central_id}: IDs are not all in parent FoF cutout")
        parent_indices = order[locations]
        central_parent = parent_indices
        details = [(central_id, len(ids))]
        remove_satellite = old_remove_satellite.copy()
    else:
        assigned = np.zeros(parent_ids.size, dtype=bool)
        for subhalo_id in candidate_subhalos(session, api_root, group_id):
            ids = subhalo_star_ids(session, api_root, subhalo_id)
            locations = np.searchsorted(sorted_parent_ids, ids)
            if np.any(locations >= sorted_parent_ids.size):
                raise RuntimeError(f"Subhalo {subhalo_id}: IDs are not all in parent FoF cutout")
            if not np.array_equal(sorted_parent_ids[locations], ids):
                raise RuntimeError(f"Subhalo {subhalo_id}: IDs are not all in parent FoF cutout")
            parent_indices = order[locations]
            if assigned[parent_indices].any():
                raise RuntimeError(f"FoF group {group_id}: overlapping subhalo particle membership")
            assigned[parent_indices] = True
            details.append((subhalo_id, len(ids)))
            if subhalo_id == central_id:
                central_parent = parent_indices
        if central_parent is None:
            raise RuntimeError(
                f"SubLink central {central_id} missing from >=100-star subhalo list "
                f"for FoF group {group_id}"
            )
        remove_satellite = assigned.copy()
        remove_satellite[central_parent] = False

    if central_parent is None:
        raise RuntimeError(f"SubLink central {central_id} has no matched star particles")
    is_central = np.zeros(parent_ids.size, dtype=bool)
    is_central[central_parent] = True
    if np.any(remove_satellite & is_central):
        raise RuntimeError(f"Central and satellite particle masks overlap in group {group_id}")
    expected_central_n = int(epoch["stellar_particle_count"])
    if int(is_central.sum()) != expected_central_n:
        raise RuntimeError(
            f"Central ParticleID count {int(is_central.sum())} differs from SubLink catalog "
            f"count {expected_central_n} for group {group_id}"
        )
    if args.snapshot == 84 and int(is_central.sum()) != old_central_nstar:
        raise RuntimeError(f"Central particle count differs from cached TNG membership for {group_id}")

    temp_fd, temp_name = tempfile.mkstemp(prefix=f".{outpath.name}.", suffix=".partial", dir=OUT)
    os.close(temp_fd)
    try:
        with h5py.File(temp_name, "w") as handle:
            group = handle.create_group("PartType4")
            for name, values in stars.items():
                group.create_dataset(name, data=values, compression="gzip", compression_opts=2,
                                     shuffle=True)
            group.create_dataset("RemoveSatelliteGalaxy", data=remove_satellite,
                                 compression="gzip", compression_opts=1)
            group.create_dataset("CentralSubhaloParticle", data=is_central,
                                 compression="gzip", compression_opts=1)
            for name, value in {
                "GroupID": group_id,
                "SampleRank": args.sample_rank,
                "Snapshot": args.snapshot,
                "Redshift": float(header["Redshift"]),
                "CentralSubhaloID": central_id,
                "NSubhalosAtLeast100Stars": len(details),
                "NParentFoFStars": int(parent_ids.size),
                "NStarsInQualifyingSubhalos": int(assigned.sum()) if args.snapshot != 84 else -1,
                "NStarsRemovedAsSatelliteGalaxies": int(remove_satellite.sum()),
                "NStarsInCentralSubhalo": int(is_central.sum()),
                "ParentCutoutBytes": parent_bytes,
                "StarFields": STAR_FIELDS,
                "MembershipMethod": (
                    ("SubLink central ParticleIDs matched exactly to the frozen 2026-10-29 "
                     "parent FoF membership file and its satellite mask" if args.snapshot == 84 else
                     "SubLink GroupFirstSub progenitor fixed as central; exact TNG "
                     "subhalo ParticleIDs matched to parent FoF PartType4 IDs; "
                     "non-central subhalos with >=100 stars removed")
                ),
            }.items():
                handle.attrs[name] = value
            subgroup = handle.create_group("QualifyingSubhalos")
            subgroup.create_dataset("SubhaloID", data=np.array([d[0] for d in details], dtype="i8"))
            subgroup.create_dataset("Nstar", data=np.array([d[1] for d in details], dtype="i8"))
            handle.flush()
        os.replace(temp_name, outpath)
    except Exception:
        try:
            os.unlink(temp_name)
        except FileNotFoundError:
            pass
        raise

    print(json.dumps({
        "sample_rank": args.sample_rank,
        "snapshot": args.snapshot,
        "group_id": group_id,
        "central_subhalo_id": central_id,
        "n_parent_type4": int(parent_ids.size),
        "n_qualifying_subhalos": len(details) if args.snapshot != 84 else -1,
        "n_central_stars": int(is_central.sum()),
        "n_satellite_stars_removed": int(remove_satellite.sum()),
        "parent_cutout_bytes": parent_bytes,
        "output": str(outpath),
        "output_bytes": outpath.stat().st_size,
    }, sort_keys=True), flush=True)


if __name__ == "__main__":
    main()
