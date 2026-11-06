#!/usr/bin/env python3
"""Render versioned growth scores with explicit isophote detection flags."""
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
    parser.add_argument("--output", type=Path,
                        default=HERE / "growth_scores_r1audit_20261106")
    args = parser.parse_args()
    if args.workers < 1 or args.workers > 4:
        raise SystemExit("Use 1-4 CPU workers")
    output_dir = args.output.resolve()
    cohort = sorted(int(row["sample_rank"]) for row in read(HERE / "lineage_cohort.csv")
                    if int(row["central_at_all_targets"]) == 1)
    output_dir.mkdir(parents=True, exist_ok=True)
    tasks = [rank for rank in cohort
             if not (output_dir / f"tree_rank_{rank:02d}_growth_scores.csv").is_file()]
    progress = output_dir.parent / f"{output_dir.name}_progress.tsv"
    with progress.open("a", newline="") as stream:
        writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
        if progress.stat().st_size == 0:
            writer.writerow(("sample_rank", "status", "detail"))
        print(json.dumps({"cohort": len(cohort), "already_rendered": len(cohort) - len(tasks),
                          "to_render": len(tasks), "workers": args.workers}), flush=True)
        with ThreadPoolExecutor(max_workers=args.workers) as executor:
            futures = {
                executor.submit(subprocess.run, [
                    sys.executable, str(RENDERER), "--sample-rank", str(rank),
                    "--output-directory", str(output_dir)],
                    capture_output=True, text=True): rank
                for rank in tasks
            }
            failures = 0
            for completed, future in enumerate(as_completed(futures), 1):
                rank = futures[future]
                result = future.result()
                output = output_dir / f"tree_rank_{rank:02d}_growth_scores.csv"
                status = "ok" if result.returncode == 0 and output.is_file() else "failed"
                if status == "failed":
                    failures += 1
                lines = result.stdout.strip().splitlines()
                detail = lines[-1] if lines else ""
                if status == "failed" and result.stderr.strip():
                    detail = result.stderr.strip().splitlines()[-1]
                writer.writerow((rank, status, detail))
                stream.flush()
                print(json.dumps({"completed": completed, "total": len(tasks),
                                  "sample_rank": rank, "status": status,
                                  "failures": failures}), flush=True)
        print(json.dumps({"finished": True, "output_directory": str(output_dir),
                          "attempted": len(tasks),
                          "failures": failures, "progress_file": str(progress)}), flush=True)
        if failures:
            raise SystemExit(1)


if __name__ == "__main__":
    main()
