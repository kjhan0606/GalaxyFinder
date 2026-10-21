#!/bin/bash
#SBATCH --job-name=darwin-newdd
#SBATCH --partition=a10
#SBATCH --nodes=2
#SBATCH --ntasks=25
#SBATCH --ntasks-per-node=13
#SBATCH --cpus-per-task=1
#SBATCH --mem=300G
#SBATCH --time=24:00:00
#SBATCH --array=0-3%1
#SBATCH --distribution=block:block
#SBATCH --chdir=/gpfs/kjhan/Others/Jaehyun/2026.Apr/Darwin1/Cocoa_analysis_20261020
#SBATCH --output=/gpfs/kjhan/Others/Jaehyun/2026.Apr/Darwin1/Cocoa_analysis_20261020/logs/newdd_%A_%a.out
#SBATCH --error=/gpfs/kjhan/Others/Jaehyun/2026.Apr/Darwin1/Cocoa_analysis_20261020/logs/newdd_%A_%a.err

set -euo pipefail
module purge
module load intel/mpi/2021.17

SNAPS=(972 977 980 984)
SNAP=${SNAPS[$SLURM_ARRAY_TASK_ID]}
printf -v PADDED '%05d' "$SNAP"

export PATH=/opt/ohpc/pub/intel/oneapi/mpi/2021.17/bin:${PATH}
export LD_LIBRARY_PATH=/opt/ohpc/pub/intel/oneapi/mpi/2021.17/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}
export I_MPI_PIN=0
# Skip AMR, gas, and sink data because this reader cannot process the output's
# AMR bookkeeping; retain DM/star particles and full stellar properties.
export NEWDD_SKIP_AMR=1
# Pass the comoving box size that NewGalFinder expects. This value was derived
# independently from unit_l, H0, and aexp in all four selected info files.
export NEWDD_BOXLEN_CMPCH=44.171267363

test -f "output_${PADDED}/info_${PADDED}.txt"
test ! -e "FoF_Data/NewDD.${PADDED}"
printf 'NewDD start snapshot=%s time=%s nodes=%s\n' "$PADDED" "$(date -Is)" "$SLURM_JOB_NODELIST"

mpirun -n 25 /gpfs/kjhan/Others/Jaehyun/2026.Apr/NewDD/newdd.exe "$SNAP" 100

test -s "FoF_Data/NewDD.${PADDED}/SN.${PADDED}.DM.00000.dat"
test -s "FoF_Data/NewDD.${PADDED}/SN.${PADDED}.STAR.00000.dat"
printf 'NewDD complete snapshot=%s time=%s\n' "$PADDED" "$(date -Is)"
