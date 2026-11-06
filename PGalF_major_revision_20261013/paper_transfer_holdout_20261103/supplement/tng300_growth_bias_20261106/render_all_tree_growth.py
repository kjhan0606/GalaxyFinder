#!/usr/bin/env python3
"""Resume-safe render and measurement pass for the 36-tree analysis sample."""
from __future__ import annotations

import argparse
import csv
import json
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

HERE = Path(__file__).resolve().parent
RENDERER = HERE / "render_tng_growth_crossrender.py"


def read(path: Path):
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--workers", type=int, default=2)
    args = parser.parse_args()
    if args.workers < 1 or args.workers > 4:
        raise SystemExit("Use 1-4 CPU workers")
    cohort = sorted(int(r["sample_rank"]) for r in read(HERE / "lineage_cohort.csv")
                    if int(r["central_at_all_targets"]) == 1)
    expected = []
    for rank in cohort:
        for snapshot in (84, 67, 50):
            item = [r for r in read(HERE / "lineage_epochs.csv")
                    if int(r["sample_rank"]) == rank and int(r["snapshot"]) == snapshot]
            if len(item) != 1:
                raise RuntimeError(f"Missing unique tree epoch for rank={rank}, snap={snapshot}")
            group_id = int(item[0]["descendant_group_id"] if snapshot == 84
                           else item[0]["progenitor_fof_group_id"])
            member = HERE / "tree_epoch_membership" / (
                f"snap{snapshot:03d}_rank{rank:02d}_group{group_id:06d}.h5")
            if not member.is_file():
                raise FileNotFoundError(f"Missing exact membership file: {member.name}")
            expected.append(member)

    OUT = HERE / "growth_scores"
    OUT.mkdir(parents=True, exist_ok=True)
    tasks = [rank for rank in cohort if not (OUT / f"tree_rank_{rank:02d}_growth_scores.csv").is_file()]
    progress = HERE / "render_all_progress.tsv"
    with progress.open("a", newline="") as stream:
        writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
        if progress.stat().st_size == 0:
            writer.writerow(("sample_rank", "status", "detail"))
        print(json.dumps({"cohort": len(cohort), "membership_files": len(expected),
                          "already_rendered": len(cohort) - len(tasks),
                          "to_render": len(tasks), "workers": args.workers}), flush=True)
        with ThreadPoolExecutor(max_workers=args.workers) as executor:
            futures = {
                executor.submit(subprocess.run, [
                    sys.executable, str(RENDERER), "--sample-rank", str(rank)
                ], capture_output=True, text=True): rank
                for rank in tasks
            }
            completed = 0
            failures = 0
            for future in as_completed(futures):
                rank = futures[future]
                result = future.result()
                completed += 1
                output = OUT / f"tree_rank_{rank:02d}_growth_scores.csv"
                status = "ok" if result.returncode == 0 and output.is_file() else "failed"
                if status == "failed":
                    failures += 1
                detail = result.stdout.strip().splitlines()[-1] if result.stdout.strip() else ""
                if status == "failed" and result.stderr.strip():
                    detail = result.stderr.strip().splitlines()[-1]
                writer.writerow((rank, status, detail))
                stream.flush()
                print(json.dumps({"completed": completed, "total": len(tasks),
                                  "sample_rank": rank, "status": status,
                                  "failures": failures}), flush=True)
        print(json.dumps({"finished": True, "attempted": len(tasks),
                          "already_rendered": len(cohort) - len(tasks),
                          "failures": failures, "progress_file": str(progress)}), flush=True)
        if failures:
            raise SystemExit(1)


if __name__ == "__main__":
    main()
