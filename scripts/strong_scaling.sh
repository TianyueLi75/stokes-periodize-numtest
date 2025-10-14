#!/bin/bash

# Set up batch job settings
#SBATCH --job-name=strong_scaling_32proc
#SBATCH --nodes=16
#SBATCH --ntasks-per-node=2
#SBATCH --cpus-per-task=32
#SBATCH --time=00:40:00
#SBATCH --mem=1000g
#SBATCH --partition=gen
#SBATCH --constraint=icelake
#SBATCH --mail-user=redacted
#SBATCH --mail-type=BEGIN,END

#SBATCH --array=1-3

WORK_DIR=~/stokes-periodize-numtest
source ${WORK_DIR}/pvfmm_modules

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export KMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export NTASK=$((${SLURM_NNODES}*${SLURM_NTASKS_PER_NODE}))
export NCORES=$((${NTASK}*${OMP_NUM_THREADS}))

cd ${WORK_DIR}

Nproc=$((SLURM_NNODES*2))
echo "Nprocess in this run: $Nproc"

# Strong scaling: tol = 1e-12, gmres_tol = 1e-8

make test2 -j &&
mpirun -n $Nproc --map-by numa:pe=${OMP_NUM_THREADS} ${WORK_DIR}/bin/test2 3 48 0 3 2000 0 1e-8 1e-12 > ${WORK_DIR}/results/2000ptcls_3peri_{$Nproc}proc_{$SLURM_ARRAY_TASK_ID}_noprecond.txt 

# mpirun -n 1 --map-by numa:pe=${OMP_NUM_THREADS} ${WORK_DIR}/bin/test2 3 48 0 3 2000 0 1e-8 1e-12 > ${WORK_DIR}/results/2000ptcls_3peri_1proc_3.txt 
# mpirun -n 2 --map-by numa:pe=${OMP_NUM_THREADS} ${WORK_DIR}/bin/test2 3 48 0 3 2000 0 1e-8 1e-12 > ${WORK_DIR}/results/2000ptcls_3peri_2proc.txt 
# mpirun -n 4 --map-by numa:pe=${OMP_NUM_THREADS} ${WORK_DIR}/bin/test2 3 48 0 3 2000 0 1e-8 1e-12 > ${WORK_DIR}/results/2000ptcls_3peri_4proc.txt 
# mpirun -n 8 --map-by numa:pe=${OMP_NUM_THREADS} ${WORK_DIR}/bin/test2 3 48 0 3 2000 0 1e-8 1e-12 > ${WORK_DIR}/results/2000ptcls_3peri_8proc.txt 
# mpirun -n 16 --map-by numa:pe=${OMP_NUM_THREADS} ${WORK_DIR}/bin/test2 3 48 0 3 2000 0 1e-8 1e-12 > ${WORK_DIR}/results/2000ptcls_3peri_16proc.txt 
# mpirun -n 32 --map-by numa:pe=${OMP_NUM_THREADS} ${WORK_DIR}/bin/test2 3 48 0 3 2000 0 1e-8 1e-12 > ${WORK_DIR}/results/2000ptcls_3peri_32proc.txt 
# mpirun -n 50 --map-by numa:pe=${OMP_NUM_THREADS} ${WORK_DIR}/bin/test2 3 48 0 3 2000 0 1e-8 1e-12 > ${WORK_DIR}/results/2000ptcls_3peri_50proc.txt 
# mpirun -n 80 --map-by numa:pe=${OMP_NUM_THREADS} ${WORK_DIR}/bin/test2 3 48 0 3 2000 0 1e-8 1e-12 > ${WORK_DIR}/results/2000ptcls_3peri_80proc.txt 
