#!/bin/bash

# Set up batch job settings
#SBATCH --job-name=test_newcode
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --time=00:30:00
#SBATCH --partition=gen
# #SBATCH --partition=ccm
#SBATCH --constraint=icelake 
#SBATCH --mail-user=tianycli@umich.edu
#SBATCH --mail-type=BEGIN,END

WORK_DIR=~/stokes-periodize-numtest
source ${WORK_DIR}/pvfmm_modules

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export KMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export NTASK=$((${SLURM_NNODES}*${SLURM_NTASKS_PER_NODE}))
export NCORES=$((${NTASK}*${OMP_NUM_THREADS}))

cd ${WORK_DIR}

# mpirun -n 1 --map-by slot:pe=${OMP_NUM_THREADS} ${WORK_DIR}/bin/test2 2 16 0 1 25 0 > ${WORK_DIR}/results/25ptcls_2_16_timing.txt
make test2_ptcl_conv -j32 &&
mpirun -n 1 --map-by numa:pe=${OMP_NUM_THREADS} ${WORK_DIR}/bin/test2_ptcl_conv 5 24 1 25 0 30000 > ${WORK_DIR}/results/25ptcls_5_24_newcode_30k_noprecond.txt

# make test_proxy_mats -B -j &&
# mpirun -n 1 --map-by slot:pe=${OMP_NUM_THREADS} ./bin/test2_ptcl_periodize_mats 30 20 > proxy_mat_30_20.txt

# make test_proxy_order -B -j &&
# mpirun -n 1 --map-by slot:pe=${OMP_NUM_THREADS} ./bin/test2_ptcl_periodize_order 4 48 1 25 0 30000 > proxy_order_30_20_debug.txt
# NOTE: above uses regular matrices, works well. 
#       below uses float128 matrices
# mpirun -n 1 --map-by slot:pe=${OMP_NUM_THREADS} ./bin/test2_ptcl_periodize_order 2 16 1 25 0 30000 > proxy_order_2_20_quadreal.txt

# make test1_dispersion -B -j &&
# mpirun -n 6 --map-by numa:pe=${OMP_NUM_THREADS} ./bin/test1_dispersion 400 64 5 > dispersion.txt

