#!/bin/bash

# Set up batch job settings
#SBATCH --job-name=strong_scaling
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --time=02:45:00
#SBATCH --partition=gen

WORK_DIR=~/stokes-periodize-numtest
source ${WORK_DIR}/pvfmm_modules

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export KMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export NTASK=$((${SLURM_NNODES}*${SLURM_NTASKS_PER_NODE}))
export NCORES=$((${NTASK}*${OMP_NUM_THREADS}))

cd ${WORK_DIR}

echo "Nprocess in this run: $NTASK"

# Strong scaling: tol = 1e-12, gmres_tol = 1e-8

make test2 -j &&
mpirun -n $NTASK --map-by slot:pe=${OMP_NUM_THREADS} ${WORK_DIR}/bin/test2 3 48 0 3 2000 0 1e-8 1e-12 > ${WORK_DIR}/results/2000ptcls_3peri_${NTASK}proc_64cpus.txt 
