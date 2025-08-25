#!/bin/bash

# Set up batch job settings
#SBATCH --job-name=weak_25ptcl
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=64
#SBATCH --time=00:25:00
#SBATCH --partition=gen

WORK_DIR=~/stokes-periodize-numtest
source ${WORK_DIR}/pvfmm_modules

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export KMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export NTASK=$((${SLURM_NNODES}*${SLURM_NTASKS_PER_NODE}))
export NCORES=$((${NTASK}*${OMP_NUM_THREADS}))

cd ${WORK_DIR}

# CHECK
# mpirun -n 1 --map-by slot:pe=${OMP_NUM_THREADS} ./bin/test2_ptcl_conv 4 32 1 25 0 30000 > 25ptcls_1peri_4_32_check.txt
# mpirun -n 1 --map-by slot:pe=${OMP_NUM_THREADS} ./bin/test2_ptcl_conv 4 32 3 25 0 30000 > 25ptcls_3peri_4_32_check.txt

make test2 -j &&
# SCRIPT FOR WEAK SCALING with 3-periodic, no accuracy check, background uniform flow
mpirun -n 1 --map-by slot:pe=${OMP_NUM_THREADS} ${WORK_DIR}/bin/test2 4 32 0 3 25 0 > ${WORK_DIR}/results/25ptcls_3peri_1nodes_2.txt 
# mpirun -n 2 --map-by slot:pe=${OMP_NUM_THREADS} ${WORK_DIR}/bin/test2 4 32 0 3 50 0 > ${WORK_DIR}/results/50ptcls_3peri_2nodes.txt 
# mpirun -n 4 --map-by slot:pe=${OMP_NUM_THREADS} ${WORK_DIR}/bin/test2 4 32 0 3 100 0 > ${WORK_DIR}/results/100ptcls_3peri_4nodes.txt 
# mpirun -n 8 --map-by slot:pe=${OMP_NUM_THREADS} ${WORK_DIR}/bin/test2 4 32 0 3 200 > ${WORK_DIR}/results/200ptcls_3peri_8nodes.txt
# mpirun -n 16 --map-by slot:pe=${OMP_NUM_THREADS} ${WORK_DIR}/bin/test2 4 32 0 3 400 > ${WORK_DIR}/results/400ptcls_3peri_16nodes.txt
