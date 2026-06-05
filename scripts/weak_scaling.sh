#!/bin/bash

# Set up batch job settings
#SBATCH --job-name=weak_scaling
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --time=00:05:00
#SBATCH --partition=gen

WORK_DIR=~/stokes-periodize-numtest
source ${WORK_DIR}/pvfmm_modules

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export KMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export NTASK=$((${SLURM_NNODES}*${SLURM_NTASKS_PER_NODE}))
export NCORES=$((${NTASK}*${OMP_NUM_THREADS}))

cd ${WORK_DIR}

# SCRIPT FOR WEAK SCALING with 3-periodic, no accuracy check, background uniform flow
# Weak scaling: tol = 1e-14, gmres_tol = 1e-9;

Nptcl=$((NTASK*25)) 
echo "Nprocess in this run: $NTASK, Nptcls in this run: $Nptcl; Number of threads: $OMP_NUM_THREADS."

make timing -j &&
mpirun -n $NTASK --report-bindings --map-by slot:pe=${OMP_NUM_THREADS} ${WORK_DIR}/bin/timing 6 64 1 3 $Nptcl 0 1e-9 1e-14 > ${WORK_DIR}/out/${Nptcl}ptcls_3peri_${NTASK}proc.txt 
