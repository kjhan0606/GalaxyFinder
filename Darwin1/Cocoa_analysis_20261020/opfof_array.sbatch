#!/bin/bash
#SBATCH --job-name=darwin-fof
#SBATCH --partition=a10
#SBATCH --nodes=3
#SBATCH --ntasks=10
#SBATCH --ntasks-per-node=4
# opFoF is MPI-only; allocate one CPU per rank.
#SBATCH --cpus-per-task=1
#SBATCH --mem=256G
#SBATCH --time=48:00:00
#SBATCH --exclude=syn02
#SBATCH --array=0-3%1
#SBATCH --distribution=block:block
#SBATCH --chdir=/gpfs/kjhan/Others/Jaehyun/2026.Apr/Darwin1/Cocoa_analysis_20261020
#SBATCH --output=/gpfs/kjhan/Others/Jaehyun/2026.Apr/Darwin1/Cocoa_analysis_20261020/logs/opfof_%A_%a.out
#SBATCH --error=/gpfs/kjhan/Others/Jaehyun/2026.Apr/Darwin1/Cocoa_analysis_20261020/logs/opfof_%A_%a.err

set -euo pipefail
module purge
module load intel/mpi/2021.17

SNAPS=(972 977 980 984)
SNAP=${SNAPS[$SLURM_ARRAY_TASK_ID]}
printf -v PADDED '%05d' "$SNAP"
GARBAGE_DIR="FoF_Garbage/Garb.${PADDED}"

# opFoF truncates BottomFaceContactHalo files on the first slab. Refuse to
# write through a symlink or into a non-empty directory from an earlier run.
if [[ -L "$GARBAGE_DIR" ]]; then
  echo "refusing symlinked FoF auxiliary-output directory: $GARBAGE_DIR" >&2
  exit 3
fi
if [[ -e "$GARBAGE_DIR" && ! -d "$GARBAGE_DIR" ]]; then
  echo "FoF auxiliary-output path exists and is not a directory: $GARBAGE_DIR" >&2
  exit 3
fi
if [[ -d "$GARBAGE_DIR" && -n "$(ls -A "$GARBAGE_DIR")" ]]; then
  echo "refusing non-empty FoF auxiliary-output directory: $GARBAGE_DIR" >&2
  exit 3
fi

export PATH=/opt/ohpc/pub/intel/oneapi/mpi/2021.17/bin:${PATH}
export LD_LIBRARY_PATH=/opt/ohpc/pub/intel/oneapi/mpi/2021.17/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}
export I_MPI_COMPATIBILITY=4
export I_MPI_PIN_DOMAIN=auto

test -s "FoF_Data/NewDD.${PADDED}/SN.${PADDED}.00000.info"
for ((i=0; i<100; i++)); do
  printf -v PART '%05d' "$i"
  test -s "FoF_Data/NewDD.${PADDED}/SN.${PADDED}.DM.${PART}.dat"
  test -e "FoF_Data/NewDD.${PADDED}/SN.${PADDED}.STAR.${PART}.dat"
done
test ! -e "FoF_Data/FoF.${PADDED}"
printf 'opFoF start snapshot=%s time=%s nodes=%s\n' "$PADDED" "$(date -Is)" "$SLURM_JOB_NODELIST"

# NewDD was run with nsplit=100; opFoF's second argument is the slab-file count.
mpirun -n 10 /gpfs/kjhan/Others/Jaehyun/2026.Apr/opFoF/opfof.exe "$SNAP" 100

test -s "FoF_Data/FoF.${PADDED}/FoF_halo_cat.${PADDED}"
test -s "FoF_Data/FoF.${PADDED}/FoF_member_particle.${PADDED}"
printf 'opFoF complete snapshot=%s time=%s\n' "$PADDED" "$(date -Is)"
