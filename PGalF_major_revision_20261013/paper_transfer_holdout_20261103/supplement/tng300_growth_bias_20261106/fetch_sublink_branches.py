#!/usr/bin/env python3
"""Fetch SubLink main-progenitor branches for the locked 40-host TNG sample.

The API key is read from ~/.config/illustris-tng/api_key. It is never written
to output, logs, or request metadata. Raw branch files and a normalized table
are written only when the destination directory does not already contain
outputs, so reruns cannot silently replace data.
"""
from __future__ import annotations

import concurrent.futures
import csv
import datetime as dt
import io
import json
import time
from pathlib import Path

import h5py
import numpy as np
import requests

HERE = Path(__file__).resolve().parent
PAPER = HERE.parents[1]
SAMPLE = PAPER / "supplement/tng300_subhalopos_r50_transfer_20261103/tng300_subhalo_pos_40.csv"
KEYFILE = Path.home() / ".config/illustris-tng/api_key"
API = "https://www.tng-project.org/api/TNG300-1/"
TARGET_REDSHIFTS = (0.19728418237600986, 0.5, 1.0)


def new_session(key: str) -> requests.Session:
    session = requests.Session()
    session.headers.update({"api-key": key, "Accept": "application/json,application/octet-stream"})
    return session


def get(session: requests.Session, url: str, *, params=None, timeout=(20, 180)):
    last_status = None
    for attempt in range(6):
        try:
            response = session.get(url, params=params, timeout=timeout)
            if response.status_code in (403, 429, 500, 502, 503, 504):
                last_status = response.status_code
                if attempt < 5:
                    time.sleep(min(2**attempt, 30))
                    continue
            if not response.ok:
                raise RuntimeError(f"TNG API returned HTTP {response.status_code}; response URL omitted")
            return response
        except requests.RequestException as exc:
            last_status = type(exc).__name__
            if attempt == 5:
                break
            time.sleep(min(2**attempt, 30))
    raise RuntimeError(f"TNG API request failed after retries ({last_status}); URL and credentials omitted")


def fetch_snapshot_redshifts(session: requests.Session) -> dict[int, float]:
    payload = get(session, API + "snapshots/").json()
    if isinstance(payload, list):
        rows = payload
    elif isinstance(payload, dict):
        rows = list(payload.get("results", []))
        url = payload.get("next")
        while url:
            page = get(session, url).json()
            if not isinstance(page, dict):
                raise RuntimeError("Unexpected paginated snapshot response")
            rows.extend(page.get("results", []))
            url = page.get("next")
    else:
        raise RuntimeError(f"Unexpected snapshot metadata type: {type(payload).__name__}")
    snapshots = {int(row["number"]): float(row["redshift"]) for row in rows}
    if len(snapshots) < 80:
        raise RuntimeError(f"Snapshot metadata is incomplete: {len(snapshots)} rows")
    return snapshots


def fetch_one(row: dict[str, str], key: str) -> tuple[dict[str, str], list[dict[str, object]], bytes]:
    rank = int(row["sample_rank"])
    subhalo_id = int(row["central_subhalo_id"])
    session = new_session(key)
    url = API + f"snapshots/84/subhalos/{subhalo_id}/sublink/mpb.hdf5"
    response = get(session, url)
    if len(response.content) > 32 * 1024 * 1024:
        raise RuntimeError(f"Unexpectedly large MPB file for rank {rank}: {len(response.content)} bytes")
    with h5py.File(io.BytesIO(response.content), "r") as tree:
        required = ("SnapNum", "SubfindID", "SubhaloID", "SubhaloGrNr", "GroupFirstSub",
                    "SubhaloMassType", "SubhaloHalfmassRadType", "SubhaloPos", "SubhaloVelDisp",
                    "SubhaloLenType", "GroupMass", "GroupPos", "Group_M_Crit200", "Group_R_Crit200")
        missing = [field for field in required if field not in tree]
        if missing:
            raise RuntimeError(f"MPB for rank {rank} lacks required fields: {missing}")
        arrays = {field: tree[field][:] for field in required}
    snaps = np.asarray(arrays["SnapNum"], dtype=int)
    sfids = np.asarray(arrays["SubfindID"], dtype=int)
    if not len(snaps) or len(np.unique(snaps)) != len(snaps):
        raise RuntimeError(f"MPB for rank {rank} has empty or duplicate snapshots")
    root = np.flatnonzero(snaps == 84)
    if len(root) != 1 or int(sfids[root[0]]) != subhalo_id:
        raise RuntimeError(f"MPB root mismatch for rank {rank}")
    # Store only catalogue information needed for target selection and later
    # sample auditing. Large branch arrays remain in the raw HDF5 file.
    rows = []
    for i, snap in enumerate(snaps):
        rows.append({
            "sample_rank": rank,
            "descendant_group_id": int(row["group_id"]),
            "descendant_subhalo_id": subhalo_id,
            "snapshot": int(snap),
            "subhalo_id": int(sfids[i]),
            "tree_subhalo_id": int(np.asarray(arrays["SubhaloID"])[i]),
            "fof_group_id": int(np.asarray(arrays["SubhaloGrNr"])[i]),
            "fof_central_subhalo_id": int(np.asarray(arrays["GroupFirstSub"])[i]),
            "is_fof_central": int(np.asarray(arrays["GroupFirstSub"])[i]) == int(sfids[i]),
            "stellar_mass_1e10_msun_h": float(np.asarray(arrays["SubhaloMassType"])[i, 4]),
            "stellar_halfmass_radius_ckpc_h": float(np.asarray(arrays["SubhaloHalfmassRadType"])[i, 4]),
            "stellar_particle_count": int(np.asarray(arrays["SubhaloLenType"])[i, 4]),
            "subhalo_velocity_dispersion_km_s": float(np.asarray(arrays["SubhaloVelDisp"])[i]),
            "fof_mass_1e10_msun_h": float(np.asarray(arrays["GroupMass"])[i]),
            "fof_mcrit200_1e10_msun_h": float(np.asarray(arrays["Group_M_Crit200"])[i]),
            "fof_rcrit200_ckpc_h": float(np.asarray(arrays["Group_R_Crit200"])[i]),
        })
    return {"sample_rank": str(rank), "subhalo_id": str(subhalo_id), "bytes": str(len(response.content))}, rows, response.content


def main() -> None:
    if not KEYFILE.is_file() or not KEYFILE.read_text().strip():
        raise RuntimeError(f"Protected TNG API key file unavailable: {KEYFILE}")
    key = KEYFILE.read_text().strip()
    with SAMPLE.open(newline="") as stream:
        sample = list(csv.DictReader(stream))
    if len(sample) != 40 or {int(row["sample_rank"]) for row in sample} != set(range(1, 41)):
        raise RuntimeError(f"Expected locked 40-host sample, found {len(sample)} rows")

    existing = list(HERE.glob("*.hdf5"))
    if existing or (HERE / "branch_nodes.csv").exists() or (HERE / "target_epochs.csv").exists():
        raise RuntimeError("Refusing to overwrite existing branch outputs")

    snapshots = fetch_snapshot_redshifts(new_session(key))
    target_snapshots = {
        f"target_{z:.6f}": min(snapshots, key=lambda snap: abs(snapshots[snap] - z))
        for z in TARGET_REDSHIFTS
    }

    fetched = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool:
        futures = [pool.submit(fetch_one, row, key) for row in sample]
        for future in concurrent.futures.as_completed(futures):
            fetched.append(future.result())
    fetched.sort(key=lambda result: int(result[0]["sample_rank"]))

    branch_rows = []
    provenance = []
    for meta, rows, content in fetched:
        rank = int(meta["sample_rank"])
        out = HERE / f"mpb_rank_{rank:02d}.hdf5"
        out.write_bytes(content)
        provenance.append({
            "sample_rank": rank,
            "descendant_subhalo_id": int(meta["subhalo_id"]),
            "endpoint": f"/snapshots/84/subhalos/{meta['subhalo_id']}/sublink/mpb.hdf5",
            "bytes": int(meta["bytes"]),
            "file": out.name,
        })
        for row in rows:
            row["redshift"] = snapshots.get(int(row["snapshot"]), float("nan"))
            branch_rows.append(row)

    with (HERE / "branch_nodes.csv").open("w", newline="") as stream:
        fields = list(branch_rows[0])
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(branch_rows)

    epochs = []
    by_rank_snap = {(int(r["sample_rank"]), int(r["snapshot"])): r for r in branch_rows}
    for label, snap in target_snapshots.items():
        for row in sample:
            rank = int(row["sample_rank"])
            node = by_rank_snap.get((rank, snap))
            epochs.append({
                "target": label,
                "target_redshift": float(label.removeprefix("target_")),
                "target_snapshot": snap,
                "snapshot_redshift": snapshots[snap],
                "sample_rank": rank,
                "descendant_group_id": int(row["group_id"]),
                "descendant_subhalo_id": int(row["central_subhalo_id"]),
                "progenitor_present": int(node is not None),
                "progenitor_subhalo_id": "" if node is None else node["subhalo_id"],
                "progenitor_fof_group_id": "" if node is None else node["fof_group_id"],
                "progenitor_is_fof_central": "" if node is None else node["is_fof_central"],
                "stellar_mass_1e10_msun_h": "" if node is None else node["stellar_mass_1e10_msun_h"],
                "stellar_halfmass_radius_ckpc_h": "" if node is None else node["stellar_halfmass_radius_ckpc_h"],
                "stellar_particle_count": "" if node is None else node["stellar_particle_count"],
            })
    with (HERE / "target_epochs.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(epochs[0]))
        writer.writeheader()
        writer.writerows(epochs)

    snapshot_rows = [{"snapshot": snap, "redshift": z} for snap, z in sorted(snapshots.items())]
    with (HERE / "snapshot_redshifts.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=["snapshot", "redshift"])
        writer.writeheader()
        writer.writerows(snapshot_rows)

    metadata = {
        "simulation": "TNG300-1",
        "tree": "SubLink main progenitor branch",
        "retrieved_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "sample": str(SAMPLE.relative_to(PAPER)),
        "n_descendants": len(sample),
        "target_redshifts": list(TARGET_REDSHIFTS),
        "target_snapshots": target_snapshots,
        "branch_files": provenance,
        "api_key_written": False,
        "signed_download_urls_written": False,
    }
    (HERE / "provenance.json").write_text(json.dumps(metadata, indent=2) + "\n")
    print(json.dumps({
        "n_descendants": len(sample),
        "branch_nodes": len(branch_rows),
        "target_snapshots": target_snapshots,
        "target_epoch_presence": {
            label: sum(int(row["progenitor_present"]) for row in epochs if row["target"] == label)
            for label in target_snapshots
        },
        "target_epoch_fof_centrals": {
            label: sum(int(row["progenitor_is_fof_central"]) for row in epochs
                       if row["target"] == label and row["progenitor_present"])
            for label in target_snapshots
        },
        "raw_branch_mib": round(sum(int(x["bytes"]) for x in provenance) / 2**20, 3),
    }, sort_keys=True))


if __name__ == "__main__":
    main()
