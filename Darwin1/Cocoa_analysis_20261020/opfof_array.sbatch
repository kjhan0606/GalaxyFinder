#!/bin/bash
#SBATCH --job-name=darwin-fof
#SBATCH --partition=a10
#SBATCH --nodes=3
#SBATCH --ntasks=10
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=6
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

export PATH=/opt/ohpc/pub/intel/oneapi/mpi/2021.17/bin:${PATH}
export LD_LIBRARY_PATH=/opt/ohpc/pub/intel/oneapi/mpi/2021.17/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}
export OMP_NUM_THREADS="$SLURM_CPUS_PER_TASK"
export OMP_PLACES=cores
export OMP_WAIT_POLICY=passive
export OMP_STACKSIZE=128M
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
