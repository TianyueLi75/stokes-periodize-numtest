#!/bin/bash

# Set up batch job settings
#SBATCH --job-name=jobname
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=16
#SBATCH --time=00:30:00 
#SBATCH --partition=gen

WORK_DIR=~/stokes-periodize-numtest
source ${WORK_DIR}/pvfmm_modules

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export KMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export NTASK=$((${SLURM_NNODES}*${SLURM_NTASKS_PER_NODE}))
export NCORES=$((${NTASK}*${OMP_NUM_THREADS}))

cd ${WORK_DIR}

# Manufactured solutions
make test_manufactured_soln -j18 &&
for n in {1..8}; do
    for m in {4..80..4}; do
        echo "n = $n, m=$m"
        mpirun -n 1 --map-by slot:pe=${OMP_NUM_THREADS} ./bin/test2_ptcl_conv $n $m 25 0 30000 >> "out/Nelem_Nf_grid.txt"
    done
done