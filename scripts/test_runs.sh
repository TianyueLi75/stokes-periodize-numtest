#!/bin/bash

# Set up batch job settings
#SBATCH --job-name=make_m032_mats
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=64
#SBATCH --time=06:15:00
#SBATCH --partition=gen

WORK_DIR=~/stokes-periodize-numtest
source ${WORK_DIR}/pvfmm_modules

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export KMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export NTASK=$((${SLURM_NNODES}*${SLURM_NTASKS_PER_NODE}))
export NCORES=$((${NTASK}*${OMP_NUM_THREADS}))

cd ${WORK_DIR}

# mpirun -n 1 --map-by slot:pe=${OMP_NUM_THREADS} ${WORK_DIR}/bin/test2_ptcl_precond_left 8 64 0 1 25 0 30000 3 > ${WORK_DIR}/results/25ptcls_8_64.txt

make test_proxy_mats -B -j &&
mpirun -n 1 --map-by slot:pe=${OMP_NUM_THREADS} ./bin/test2_ptcl_periodize_mats 30 32 > proxy_mat_30_32.txt

