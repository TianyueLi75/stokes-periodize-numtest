#!/bin/bash

# Set up batch job settings
#SBATCH --job-name=weak_10proc
#SBATCH --nodes=5
#SBATCH --ntasks-per-node=2
#SBATCH --cpus-per-task=32
#SBATCH --mem=200g
#SBATCH --time=00:25:00
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

make test2 -j &&
# SCRIPT FOR WEAK SCALING with 3-periodic, no accuracy check, background uniform flow
# Weak scaling: tol = 1e-14, gmres_tol = 1e-9; New code with 60k particles and more detailed param grid gives: 3-24 is enough for 1e-6 error

# mpirun -n 1 --map-by numa:pe=${OMP_NUM_THREADS} ${WORK_DIR}/bin/test2 3 24 1 3 25 0 1e-9 1e-14 > ${WORK_DIR}/results/25ptcls_3peri_1proc_streamline.txt 
# mpirun -n 2 --map-by numa:pe=${OMP_NUM_THREADS} ${WORK_DIR}/bin/test2 3 24 0 3 50 0 1e-9 1e-14 > ${WORK_DIR}/results/50ptcls_3peri_2proc.txt 
# mpirun -n 4 --map-by numa:pe=${OMP_NUM_THREADS} ${WORK_DIR}/bin/test2 3 24 0 3 100 0 1e-9 1e-14 > ${WORK_DIR}/results/100ptcls_3peri_4proc.txt 
# mpirun -n 8 --map-by numa:pe=${OMP_NUM_THREADS} ${WORK_DIR}/bin/test2 3 24 0 3 200 0 1e-9 1e-14 > ${WORK_DIR}/results/200ptcls_3peri_8proc.txt
mpirun -n 10 --map-by numa:pe=${OMP_NUM_THREADS} ${WORK_DIR}/bin/test2 3 24 0 3 250 0 1e-9 1e-14 > ${WORK_DIR}/results/250ptcls_3peri_10proc.txt 
# mpirun -n 14 --map-by numa:pe=${OMP_NUM_THREADS} ${WORK_DIR}/bin/test2 3 24 0 3 350 0 1e-9 1e-14 > ${WORK_DIR}/results/350ptcls_3peri_14proc.txt 
# mpirun -n 16 --map-by numa:pe=${OMP_NUM_THREADS} ${WORK_DIR}/bin/test2 3 24 0 3 400 0 1e-9 1e-14 > ${WORK_DIR}/results/400ptcls_3peri_16proc.txt
# mpirun -n 20 --map-by numa:pe=${OMP_NUM_THREADS} ${WORK_DIR}/bin/test2 3 24 0 3 500 0 1e-9 1e-14 > ${WORK_DIR}/results/500ptcls_3peri_20proc.txt  
# mpirun -n 24 --map-by numa:pe=${OMP_NUM_THREADS} ${WORK_DIR}/bin/test2 3 24 0 3 600 0 1e-9 1e-14 > ${WORK_DIR}/results/600ptcls_3peri_24proc.txt 
# mpirun -n 28 --map-by numa:pe=${OMP_NUM_THREADS} ${WORK_DIR}/bin/test2 3 24 0 3 700 0 1e-9 1e-14 > ${WORK_DIR}/results/700ptcls_3peri_28proc.txt 
# mpirun -n 32 --map-by numa:pe=${OMP_NUM_THREADS} ${WORK_DIR}/bin/test2 3 24 0 3 800 0 1e-9 1e-14 > ${WORK_DIR}/results/800ptcls_3peri_32proc.txt 
# mpirun -n 40 --map-by numa:pe=${OMP_NUM_THREADS} ${WORK_DIR}/bin/test2 3 24 0 3 1000 0 1e-9 1e-14 > ${WORK_DIR}/results/1000ptcls_3peri_40proc.txt 
# mpirun -n 50 --map-by numa:pe=${OMP_NUM_THREADS} ${WORK_DIR}/bin/test2 3 24 0 3 1250 0 1e-9 1e-14 > ${WORK_DIR}/results/1250ptcls_3peri_50proc.txt 
# mpirun -n 64 --map-by numa:pe=${OMP_NUM_THREADS} ${WORK_DIR}/bin/test2 3 24 0 3 1600 0 1e-9 1e-14 > ${WORK_DIR}/results/1600ptcls_3peri_64proc.txt 
# mpirun -n 70 --map-by numa:pe=${OMP_NUM_THREADS} ${WORK_DIR}/bin/test2 3 24 0 3 1750 0 1e-9 1e-14 > ${WORK_DIR}/results/1750ptcls_3peri_70proc.txt 
# mpirun -n 76 --map-by numa:pe=${OMP_NUM_THREADS} ${WORK_DIR}/bin/test2 3 24 0 3 1900 0 1e-9 1e-14 > ${WORK_DIR}/results/1900ptcls_3peri_76proc.txt 
# mpirun -n 80 --map-by numa:pe=${OMP_NUM_THREADS} ${WORK_DIR}/bin/test2 3 24 0 3 2000 0 1e-9 1e-14 > ${WORK_DIR}/results/2000ptcls_3peri_80proc.txt 
