#!/bin/bash

# Set up batch job settings
#SBATCH --job-name=weak_2000ptcl_streamline
#SBATCH --nodes=10
#SBATCH --ntasks-per-node=2
#SBATCH --cpus-per-task=32
#SBATCH --time=03:45:00
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
# Weak scaling: tol = 1e-14, gmres_tol = 1e-9

# mpirun -n 1 --map-by numa:pe=${OMP_NUM_THREADS} ${WORK_DIR}/bin/test2 4 32 1 3 25 0 1e-9 1e-14 > ${WORK_DIR}/results/25ptcls_3peri_1nodes_numa_streamline.txt 
# mpirun -n 2 --map-by numa:pe=${OMP_NUM_THREADS} ${WORK_DIR}/bin/test2 4 32 0 3 50 0 1e-9 1e-14 > ${WORK_DIR}/results/50ptcls_3peri_2nodes_numa.txt 
# mpirun -n 4 --map-by numa:pe=${OMP_NUM_THREADS} ${WORK_DIR}/bin/test2 4 32 0 3 100 0 1e-9 1e-14 > ${WORK_DIR}/results/100ptcls_3peri_4nodes_numa.txt 
# mpirun -n 8 --map-by numa:pe=${OMP_NUM_THREADS} ${WORK_DIR}/bin/test2 4 32 0 3 200 0 1e-9 1e-14 > ${WORK_DIR}/results/200ptcls_3peri_8nodes_numa_2.txt
# mpirun -n 8 --map-by numa:pe=${OMP_NUM_THREADS} ${WORK_DIR}/bin/test2 4 32 1 3 400 0 1e-9 1e-14 > ${WORK_DIR}/results/400ptcls_3peri_8nodes_numa_streamline.txt 
# mpirun -n 32 --map-by numa:pe=${OMP_NUM_THREADS} ${WORK_DIR}/bin/test2 4 32 0 3 800 0 1e-9 1e-14 > ${WORK_DIR}/results/800ptcls_3peri_32nodes_numa.txt 
# mpirun -n 64 --map-by numa:pe=${OMP_NUM_THREADS} ${WORK_DIR}/bin/test2 4 32 0 3 1600 0 1e-9 1e-14 > ${WORK_DIR}/results/1600ptcls_3peri_64nodes_numa.txt 
# mpirun -n 20 --map-by numa:pe=${OMP_NUM_THREADS} ${WORK_DIR}/bin/test2 4 32 1 3 2000 0 1e-9 1e-14 > ${WORK_DIR}/results/2000ptcls_3peri_20nodes_numa.txt 

# Gridded spheres starting with 27ptcls and up
# mpirun -n 1 --map-by numa:pe=${OMP_NUM_THREADS} ${WORK_DIR}/bin/test2 4 32 0 3 27 0 1e-9 1e-14 > ${WORK_DIR}/results/27ptcls_3peri_1nodes.txt 
# mpirun -n 4 --map-by numa:pe=${OMP_NUM_THREADS} ${WORK_DIR}/bin/test2 4 32 0 3 64 0 1e-9 1e-14 > ${WORK_DIR}/results/64ptcls_3peri_4nodes.txt 
# mpirun -n 6 --map-by numa:pe=${OMP_NUM_THREADS} ${WORK_DIR}/bin/test2 4 32 0 3 125 0 1e-9 1e-14 > ${WORK_DIR}/results/125ptcls_3peri_6nodes.txt 
# mpirun -n 14 --map-by numa:pe=${OMP_NUM_THREADS} ${WORK_DIR}/bin/test2 4 32 0 3 343 0 1e-9 1e-14 > ${WORK_DIR}/results/343ptcls_3peri_14nodes.txt
# mpirun -n 20 --map-by numa:pe=${OMP_NUM_THREADS} ${WORK_DIR}/bin/test2 4 32 0 3 512 0 1e-9 1e-14 > ${WORK_DIR}/results/512ptcls_3peri_20nodes.txt
# mpirun -n 28 --map-by numa:pe=${OMP_NUM_THREADS} ${WORK_DIR}/bin/test2 4 32 0 3 729 0 1e-9 1e-14 > ${WORK_DIR}/results/729ptcls_3peri_28nodes.txt 
# mpirun -n 38 --map-by numa:pe=${OMP_NUM_THREADS} ${WORK_DIR}/bin/test2 4 32 0 3 1000 0 1e-9 1e-14 > ${WORK_DIR}/results/1000ptcls_3peri_38nodes.txt
# mpirun -n 50 --map-by numa:pe=${OMP_NUM_THREADS} ${WORK_DIR}/bin/test2 4 32 0 3 1331 0 1e-9 1e-14 > ${WORK_DIR}/results/1331ptcls_3peri_50nodes.txt
# mpirun -n 64 --map-by numa:pe=${OMP_NUM_THREADS} ${WORK_DIR}/bin/test2 4 32 0 3 1728 0 1e-9 1e-14 > ${WORK_DIR}/results/1728ptcls_3peri_64nodes.txt
mpirun -n 80 --map-by numa:pe=${OMP_NUM_THREADS} ${WORK_DIR}/bin/test2 4 32 0 3 2197 0 1e-9 1e-14 > ${WORK_DIR}/results/2197ptcls_3peri_80nodes_2.txt 
