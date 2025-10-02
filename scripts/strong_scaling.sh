#!/bin/bash

# Set up batch job settings
#SBATCH --job-name=strong_scaling_test1nodes
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --time=03:40:00
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

# Strong scaling: tol = 1e-12, gmres_tol = 1e-8

make test2 -j &&
mpirun -n 1 --map-by numa:pe=${OMP_NUM_THREADS} ${WORK_DIR}/bin/test2 2 20 0 3 800 0 1e-8 1e-12 > ${WORK_DIR}/results/800ptcls_3peri_1nodes_numa_2.txt 
# mpirun -n 2 --map-by numa:pe=${OMP_NUM_THREADS} ${WORK_DIR}/bin/test2 2 20 0 3 800 0 1e-8 1e-12 > ${WORK_DIR}/results/800ptcls_3peri_2nodes_numa_2.txt 
# mpirun -n 4 --map-by numa:pe=${OMP_NUM_THREADS} ${WORK_DIR}/bin/test2 2 20 0 3 800 0 1e-8 1e-12 > ${WORK_DIR}/results/800ptcls_3peri_4nodes_numa_2.txt 
# mpirun -n 8 --map-by numa:pe=${OMP_NUM_THREADS} ${WORK_DIR}/bin/test2 2 20 0 3 800 0 1e-8 1e-12 > ${WORK_DIR}/results/800ptcls_3peri_8nodes_numa_2.txt 
# mpirun -n 16 --map-by numa:pe=${OMP_NUM_THREADS} ${WORK_DIR}/bin/test2 2 20 0 3 800 0 1e-8 1e-12 > ${WORK_DIR}/results/800ptcls_3peri_16nodes_numa_2.txt 
# mpirun -n 32 --map-by numa:pe=${OMP_NUM_THREADS} ${WORK_DIR}/bin/test2 2 20 0 3 800 0 1e-8 1e-12 > ${WORK_DIR}/results/800ptcls_3peri_32nodes_numa_2.txt 
# mpirun -n 50 --map-by numa:pe=${OMP_NUM_THREADS} ${WORK_DIR}/bin/test2 2 20 0 3 800 0 1e-8 1e-12 > ${WORK_DIR}/results/800ptcls_3peri_50nodes_numa_2.txt 
# mpirun -n 80 --map-by numa:pe=${OMP_NUM_THREADS} ${WORK_DIR}/bin/test2 2 20 0 3 800 0 1e-8 1e-12 > ${WORK_DIR}/results/800ptcls_3peri_80nodes_numa_2.txt 
