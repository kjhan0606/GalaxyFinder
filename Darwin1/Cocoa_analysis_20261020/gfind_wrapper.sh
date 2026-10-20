#!/bin/bash
set -euo pipefail

RANK="${SLURM_PROCID:-${PMI_RANK:-${PMIX_RANK:-${OMPI_COMM_WORLD_RANK:-0}}}}"
case "$RANK" in
  0) export OMP_NUM_THREADS=1; CORES=0 ;;
  1) export OMP_NUM_THREADS=20; CORES=1-20 ;;
  2) export OMP_NUM_THREADS=20; CORES=21-40 ;;
  3) export OMP_NUM_THREADS=20; CORES=1-20 ;;
  4) export OMP_NUM_THREADS=20; CORES=21-40 ;;
  *) echo "unexpected MPI rank $RANK" >&2; exit 2 ;;
esac
export OMP_PROC_BIND=close
export OMP_PLACES=cores
printf 'rank=%s host=%s omp=%s cores=%s snapshot=%s\n' "$RANK" "$(hostname)" "$OMP_NUM_THREADS" "$CORES" "$GFINDSNAP"
exec taskset -c "$CORES" /gpfs/kjhan/Others/Jaehyun/2026.Apr/NewGalFinder/gfind.exe "$GFINDSNAP"
