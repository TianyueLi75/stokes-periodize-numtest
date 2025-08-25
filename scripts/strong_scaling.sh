#!/bin/bash

# Set up batch job settings
#SBATCH --job-name=strong_scaling_3peri_2nodes
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=2
#SBATCH --cpus-per-task=15
#SBATCH --time=01:00:00
#SBATCH --partition=gen

WORK_DIR=~/stokes-periodize-numtest
source ${WORK_DIR}/pvfmm_modules

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export KMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export NTASK=$((${SLURM_NNODES}*${SLURM_NTASKS_PER_NODE}))
export NCORES=$((${NTASK}*${OMP_NUM_THREADS}))

cd ${WORK_DIR}

# SCRIPT FOR 400 PTCLS STRONG SCALING (fix problem, increase processes)
# mpirun --report-bindings -n 1 --map-by slot:pe=${OMP_NUM_THREADS} ${WORK_DIR}/bin/test2_ptcl_precond_left 2 32 0 3 400 0 30000 3 > ${WORK_DIR}/results/400ptcls_3peri_1nodes_1K1A.txt 
mpirun --report-bindings -n 2 --map-by slot:pe=${OMP_NUM_THREADS} ${WORK_DIR}/bin/test2_ptcl_precond_left 2 32 0 3 400 0 30000 3 > ${WORK_DIR}/results/400ptcls_3peri_2nodes_1K1A_2.txt 
# mpirun --report-bindings -n 4 --map-by slot:pe=${OMP_NUM_THREADS} ${WORK_DIR}/bin/test2_ptcl_precond_left 2 32 0 3 400 0 30000 3 > ${WORK_DIR}/results/400ptcls_3peri_4nodes_1K1A.txt 
# mpirun --report-bindings -n 8 --map-by slot:pe=${OMP_NUM_THREADS} ${WORK_DIR}/bin/test2_ptcl_precond_left 2 32 0 3 400 0 30000 3 > ${WORK_DIR}/results/400ptcls_3peri_8nodes_1K1A.txt 
