#!/bin/bash
#SBATCH --job-name=darwin-gal972
#SBATCH --partition=a100
#SBATCH --nodes=2
#SBATCH --ntasks=5
#SBATCH --ntasks-per-node=3
#SBATCH --cpus-per-task=20
#SBATCH --gres=gpu:1
#SBATCH --mem=600G
#SBATCH --time=48:00:00
#SBATCH --distribution=block:block
#SBATCH --chdir=/gpfs/kjhan/Others/Jaehyun/2026.Apr/Darwin1/Cocoa_analysis_20261020
#SBATCH --output=/gpfs/kjhan/Others/Jaehyun/2026.Apr/Darwin1/Cocoa_analysis_20261020/logs/gfind_00972_%j.out
#SBATCH --error=/gpfs/kjhan/Others/Jaehyun/2026.Apr/Darwin1/Cocoa_analysis_20261020/logs/gfind_00972_%j.err

set -euo pipefail
module purge
module load intel/mpi/2021.17

export PATH=/opt/ohpc/pub/intel/oneapi/mpi/2021.17/bin:${PATH}
export LD_LIBRARY_PATH=/opt/ohpc/pub/intel/oneapi/mpi/2021.17/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}
export I_MPI_PIN=0
export GFINDSNAP=972

test -s FoF_Data/FoF.00972/FoF_halo_cat.00972
test -s FoF_Data/FoF.00972/FoF_member_particle.00972
test ! -e FoF_Data/FoF.00972/GALCATALOG.LIST.00972
printf 'galaxy finder start snapshot=00972 time=%s nodes=%s\n' "$(date -Is)" "$SLURM_JOB_NODELIST"

mpirun -n 5 -ppn 3 ./gfind_wrapper.sh

test -s FoF_Data/FoF.00972/GALCATALOG.LIST.00972
test -s FoF_Data/FoF.00972/GALFIND.DATA.00972
printf 'galaxy finder complete snapshot=00972 time=%s\n' "$(date -Is)"
