#!/bin/bash
#SBATCH --job-name=darwin-fof972
#SBATCH --partition=a10
#SBATCH --nodes=3
#SBATCH --ntasks=10
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=6
#SBATCH --gres=gpu:1
#SBATCH --mem=256G
#SBATCH --time=48:00:00
#SBATCH --exclude=syn02
#SBATCH --distribution=block:block
#SBATCH --chdir=/gpfs/kjhan/Others/Jaehyun/2026.Apr/Darwin1/Cocoa_analysis_20261020
#SBATCH --output=/gpfs/kjhan/Others/Jaehyun/2026.Apr/Darwin1/Cocoa_analysis_20261020/logs/opfof_00972_%j.out
#SBATCH --error=/gpfs/kjhan/Others/Jaehyun/2026.Apr/Darwin1/Cocoa_analysis_20261020/logs/opfof_00972_%j.err

set -euo pipefail
module purge
module load intel/mpi/2021.17

export PATH=/opt/ohpc/pub/intel/oneapi/mpi/2021.17/bin:${PATH}
export LD_LIBRARY_PATH=/opt/ohpc/pub/intel/oneapi/mpi/2021.17/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}
export OMP_NUM_THREADS="$SLURM_CPUS_PER_TASK"
export OMP_PLACES=cores
export OMP_WAIT_POLICY=passive
export OMP_STACKSIZE=128M
export I_MPI_COMPATIBILITY=4
export I_MPI_PIN_DOMAIN=auto

test -s FoF_Data/NewDD.00972/SN.00972.DM.00000.dat
test ! -e FoF_Data/FoF.00972
for i in 00000 00001 00002 00003 00004 00005 00006 00007 00008 00009; do
  test -s "FoF_Garbage/Garb.00972/BottomFaceContactHalo.${i}.00972.dat"
done
printf 'opFoF start snapshot=00972 time=%s nodes=%s\n' "$(date -Is)" "$SLURM_JOB_NODELIST"

mpirun -n 10 /gpfs/kjhan/Others/Jaehyun/2026.Apr/opFoF/opfof.exe 972 200

test -s FoF_Data/FoF.00972/FoF_halo_cat.00972
test -s FoF_Data/FoF.00972/FoF_member_particle.00972
printf 'opFoF complete snapshot=00972 time=%s\n' "$(date -Is)"
