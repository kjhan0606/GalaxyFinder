#!/usr/bin/env python3
"""Resume-safe concurrent fetch of the locked 36-tree, three-epoch sample."""
from __future__ import annotations

import argparse
import csv
import json
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

HERE = Path(__file__).resolve().parent
FETCHER = HERE / "fetch_tree_epoch_membership.py"


def rows(path: Path):
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--workers", type=int, default=3)
    parser.add_argument("--snapshot84-parent-dir", type=Path, default=None,
                        help="Directory containing original snapshot-84 FoF membership cutouts")
    args = parser.parse_args()
    if args.workers < 1 or args.workers > 4:
        raise SystemExit("Use 1-4 concurrent API workers")
    cohort = [r for r in rows(HERE / "lineage_cohort.csv")
              if int(r["central_at_all_targets"]) == 1]
    epochs = rows(HERE / "lineage_epochs.csv")
    tasks = []
    skipped = 0
    for snap in (84, 67, 50):
        for tree in cohort:
            rank = int(tree["sample_rank"])
            epoch = [r for r in epochs if int(r["sample_rank"]) == rank
                     and int(r["snapshot"]) == snap]
            if len(epoch) != 1:
                raise RuntimeError(f"Missing unique epoch row for rank={rank}, snapshot={snap}")
            gid = int(epoch[0]["progenitor_fof_group_id"])
            output = HERE / "tree_epoch_membership" / (
                f"snap{snap:03d}_rank{rank:02d}_group{gid:06d}.h5")
            if output.is_file():
                skipped += 1
            else:
                tasks.append((snap, rank, gid, output))

    progress = HERE / "fetch_all_progress.tsv"
    with progress.open("a", newline="") as stream:
        writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
        if progress.stat().st_size == 0:
            writer.writerow(("snapshot", "sample_rank", "group_id", "status", "detail"))
        print(json.dumps({"cohort": len(cohort), "target_files": len(cohort) * 3,
                          "already_present": skipped, "to_fetch": len(tasks),
                          "workers": args.workers}), flush=True)
        with ThreadPoolExecutor(max_workers=args.workers) as executor:
            futures = {
                executor.submit(subprocess.run, [
                    sys.executable, str(FETCHER), "--snapshot", str(snap),
                    "--sample-rank", str(rank),
                    *( ["--snapshot84-parent-dir", str(args.snapshot84_parent_dir.resolve())]
                       if snap == 84 and args.snapshot84_parent_dir else [] )
                ], capture_output=True, text=True): (snap, rank, gid, output)
                for snap, rank, gid, output in tasks
            }
            completed = 0
            failures = 0
            for future in as_completed(futures):
                snap, rank, gid, output = futures[future]
                result = future.result()
                completed += 1
                status = "ok" if result.returncode == 0 and output.is_file() else "failed"
                if status == "failed":
                    failures += 1
                detail = result.stdout.strip().splitlines()[-1] if result.stdout.strip() else ""
                if status == "failed" and result.stderr.strip():
                    # The fetcher sanitizes request errors; retain only the final diagnostic line.
                    detail = result.stderr.strip().splitlines()[-1]
                writer.writerow((snap, rank, gid, status, detail))
                stream.flush()
                print(json.dumps({"completed": completed, "total": len(tasks),
                                  "snapshot": snap, "sample_rank": rank,
                                  "status": status, "failures": failures}), flush=True)
        print(json.dumps({"finished": True, "attempted": len(tasks),
                          "already_present": skipped, "failures": failures,
                          "progress_file": str(progress)}), flush=True)
        if failures:
            raise SystemExit(1)


if __name__ == "__main__":
    main()
