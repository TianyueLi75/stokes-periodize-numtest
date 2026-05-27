#!/bin/bash

# Set up batch job settings
#SBATCH --job-name=25ptcl_selfconv
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --time=01:35:00
# #SBATCH --mem=200g
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

# make vslip -j &&
# mpirun -n 1 --map-by slot:pe=${OMP_NUM_THREADS} ${WORK_DIR}/bin/test_vslip 1 1

# test periodicity
# make test_peri -j &&
# mpirun -n 1 --map-by slot:pe=${OMP_NUM_THREADS} ${WORK_DIR}/bin/test_periodicity 2 32 1 1 1e-8 1e-8 > ${WORK_DIR}/out/1ptcls_1peri_periodicity.txt

# make test_selfconv -j &&
# mpirun -n 1 --map-by slot:pe=${OMP_NUM_THREADS} ${WORK_DIR}/bin/test_convergence 1 3 25 0 > ${WORK_DIR}/out/test_convergence_25ptcls_3peri.txt
# mpirun -n ${NTASK} --map-by slot:pe=${OMP_NUM_THREADS} ${WORK_DIR}/bin/test_convergence 2 1 > ${WORK_DIR}/out/test_convergence_convdiv_1.txt 
# mpirun -n 1 --map-by slot:pe=${OMP_NUM_THREADS} ${WORK_DIR}/bin/test_convergence 4 2 3 0 > ${WORK_DIR}/out/test_convergence_plane_3spheres.txt
# mpirun -n ${NTASK} --map-by slot:pe=${OMP_NUM_THREADS} ${WORK_DIR}/bin/test_convergence 0 1 > ${WORK_DIR}/out/test_convergence_trefoil.txt

# make examples -j &&
# mpirun -n ${NTASK} --map-by slot:pe=${OMP_NUM_THREADS} ${WORK_DIR}/bin/examples 0 > ${WORK_DIR}/out/channel_example_noslip.txt
# mpirun -n ${NTASK} --map-by slot:pe=${OMP_NUM_THREADS} ${WORK_DIR}/bin/examples 4 > ${WORK_DIR}/out/plane_example_3.txt
# mpirun -n 1 --map-by slot:pe=${OMP_NUM_THREADS} ${WORK_DIR}/bin/examples 2 > ${WORK_DIR}/out/particle_example.txt

# make dispersion -j &&
# mpirun -n ${NTASK} --map-by slot:pe=${OMP_NUM_THREADS} ${WORK_DIR}/bin/dispersion 400 64 > ${WORK_DIR}/out/trefoil_vis.txt

# Streamlines
# make timing -j32 &&
# mpirun -n 1 --map-by numa:pe=${OMP_NUM_THREADS} ${WORK_DIR}/bin/test2 4 32 1 3 25 0 1e-9 1e-14 > ${WORK_DIR}/results/25ptcls_3peri_1proc_streamlines.txt
# mpirun -n 16 --map-by numa:pe=${OMP_NUM_THREADS} ${WORK_DIR}/bin/test2 4 32 1 3 400 0 1e-9 1e-14 > ${WORK_DIR}/results/400ptcls_3peri_16proc_streamlines.txt
# mpirun -n 80 --map-by numa:pe=${OMP_NUM_THREADS} ${WORK_DIR}/bin/test2 4 32 1 3 2000 0 1e-9 1e-14 > ${WORK_DIR}/results/2000ptcls_3peri_80proc_streamlines.txt

