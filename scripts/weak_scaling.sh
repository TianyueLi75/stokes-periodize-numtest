#!/bin/bash

# Set up batch job settings
#SBATCH --job-name=weak_scaling
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=100g
#SBATCH --time=00:05:00
# #SBATCH --partition=gen
#SBATCH --partition=gen
#SBATCH --constraint=icelake
#SBATCH --mail-user=redacted
#SBATCH --mail-type=BEGIN,END

WORK_DIR=~/stokes-periodize-numtest
source ${WORK_DIR}/pvfmm_modules

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export KMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export NTASK=$((${SLURM_NNODES}*${SLURM_NTASKS_PER_NODE}))
export NCORES=$((${NTASK}*${OMP_NUM_THREADS}))

cd ${WORK_DIR}

# SCRIPT FOR WEAK SCALING with 3-periodic, no accuracy check, background uniform flow
# Weak scaling: tol = 1e-14, gmres_tol = 1e-9; New code with 60k particles and more detailed param grid gives: 3-24 is enough for 1e-6 error

# use "grep -n ^+-S 50ptcls_3peri_2proc.txt" to get line number with top profiling results. Can then pipe to grab first value if needed

Nptcl=$((NTASK*25)) 
# tODO: 1) change this for when in-between Nptcls. 2) Will only go up to 32 processes = 800 particles. Request larger job?
echo "Nprocess in this run: $NTASK, Nptcls in this run: $Nptcl; Number of threads: $OMP_NUM_THREADS."

make test2 -j &&
mpirun -n $NTASK --report-bindings --map-by slot:pe=${OMP_NUM_THREADS} ${WORK_DIR}/bin/test2 2 32 1 3 $Nptcl 0 1e-6 1e-6 > ${WORK_DIR}/results/${Nptcl}ptcls_3peri_${NTASK}proc_noslip_testcode.txt 
