#!/bin/bash
#SBATCH --job-name=darwin-gal
#SBATCH --partition=a100
#SBATCH --nodes=2
#SBATCH --ntasks=5
#SBATCH --ntasks-per-node=3
#SBATCH --cpus-per-task=20
#SBATCH --mem=600G
#SBATCH --time=48:00:00
#SBATCH --array=0-3%1
#SBATCH --distribution=block:block
#SBATCH --chdir=/gpfs/kjhan/Others/Jaehyun/2026.Apr/Darwin1/Cocoa_analysis_20261020
#SBATCH --output=/gpfs/kjhan/Others/Jaehyun/2026.Apr/Darwin1/Cocoa_analysis_20261020/logs/gfind_%A_%a.out
#SBATCH --error=/gpfs/kjhan/Others/Jaehyun/2026.Apr/Darwin1/Cocoa_analysis_20261020/logs/gfind_%A_%a.err

set -euo pipefail
module purge
module load intel/mpi/2021.17

SNAPS=(972 977 980 984)
SNAP=${SNAPS[$SLURM_ARRAY_TASK_ID]}
printf -v PADDED '%05d' "$SNAP"

export PATH=/opt/ohpc/pub/intel/oneapi/mpi/2021.17/bin:${PATH}
export LD_LIBRARY_PATH=/opt/ohpc/pub/intel/oneapi/mpi/2021.17/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}
export I_MPI_PIN=0
export GFINDSNAP="$SNAP"

test -s "FoF_Data/FoF.${PADDED}/FoF_halo_cat.${PADDED}"
test -s "FoF_Data/FoF.${PADDED}/FoF_member_particle.${PADDED}"
test ! -e "FoF_Data/FoF.${PADDED}/GALCATALOG.LIST.${PADDED}"
printf 'galaxy finder start snapshot=%s time=%s nodes=%s\n' "$PADDED" "$(date -Is)" "$SLURM_JOB_NODELIST"

mpirun -n 5 -ppn 3 ./gfind_wrapper.sh

test -s "FoF_Data/FoF.${PADDED}/GALCATALOG.LIST.${PADDED}"
test -s "FoF_Data/FoF.${PADDED}/GALFIND.DATA.${PADDED}"
printf 'galaxy finder complete snapshot=%s time=%s\n' "$PADDED" "$(date -Is)"
