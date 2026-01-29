#!/bin/bash

# Set up batch job settings
#SBATCH --job-name=plane_selfconv
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=2
#SBATCH --cpus-per-task=32
#SBATCH --time=01:30:00
# #SBATCH --mem=50g
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

# test periodicity
# make test_plane -j &&
# mpirun -n 1 --map-by numa:pe=${OMP_NUM_THREADS} ${WORK_DIR}/bin/test_plane 40 2 0.01 0.9 1e-8 1e-8 > ${WORK_DIR}/out/test_plane_regluarized.txt

# make test_peri -j &&
# mpirun -n 1 --map-by slot:pe=${OMP_NUM_THREADS} ${WORK_DIR}/bin/test_periodicity 2 32 1 1 1e-8 1e-8 > ${WORK_DIR}/out/1ptcls_1peri_periodicity.txt

# make test3 -j &&
# mpirun -n 4 --map-by slot:pe=${OMP_NUM_THREADS} ${WORK_DIR}/bin/test3 > ${WORK_DIR}/out/channel_example.txt

make test_selfconv -B -j &&
# mpirun -n ${NTASK} --map-by slot:pe=${OMP_NUM_THREADS} ${WORK_DIR}/bin/test_convergence 1 2 > ${WORK_DIR}/out/test_convergence_3ptcls_2peri.txt
# mpirun -n ${NTASK} --map-by slot:pe=${OMP_NUM_THREADS} ${WORK_DIR}/bin/test_convergence 1 3 > ${WORK_DIR}/out/test_convergence_3ptcls_3peri.txt
mpirun -n ${NTASK} --map-by slot:pe=${OMP_NUM_THREADS} ${WORK_DIR}/bin/test_convergence 2 1 > ${WORK_DIR}/out/test_convergence_convdiv.txt 
mpirun -n 1 --map-by slot:pe=${OMP_NUM_THREADS} ${WORK_DIR}/bin/test_convergence 4 2 > ${WORK_DIR}/out/test_convergence_plane_3loops.txt
# mpirun -n ${NTASK} --map-by slot:pe=${OMP_NUM_THREADS} ${WORK_DIR}/bin/test_convergence 0 1 > ${WORK_DIR}/out/test_convergence_trefoil.txt

# make examples -j &&
# mpirun -n ${NTASK} --map-by slot:pe=${OMP_NUM_THREADS} ${WORK_DIR}/bin/dispersion 400 64 > ${WORK_DIR}/out/channel_example.txt

# make dispersion -j &&
# mpirun -n ${NTASK} --map-by slot:pe=${OMP_NUM_THREADS} ${WORK_DIR}/bin/dispersion 400 64 > ${WORK_DIR}/out/trefoil_vis.txt

# Timing for 25 particles
# mpirun -n 1 --map-by slot:pe=${OMP_NUM_THREADS} ${WORK_DIR}/bin/test2_timing 2 16 0 1 25 0 > ${WORK_DIR}/results/25ptcls_2_16_timing.txt
# make test2_ptcl_conv -j32 &&
# mpirun -n 1 --map-by numa:pe=${OMP_NUM_THREADS} ${WORK_DIR}/bin/test2_ptcl_conv 5 24 25 0 30000 > ${WORK_DIR}/results/25ptcls_5_24_newcode_30k.txt

# Streamlines
# make test2_timing -j32 &&
# mpirun -n 1 --map-by numa:pe=${OMP_NUM_THREADS} ${WORK_DIR}/bin/test2_timing 3 24 1 3 25 0 1e-9 1e-14 > ${WORK_DIR}/results/25ptcls_3peri_1proc_streamlines.txt
# mpirun -n 16 --map-by numa:pe=${OMP_NUM_THREADS} ${WORK_DIR}/bin/test2_timing 3 24 1 3 400 0 1e-9 1e-14 > ${WORK_DIR}/results/400ptcls_3peri_16proc_streamlines.txt
# mpirun -n 80 --map-by numa:pe=${OMP_NUM_THREADS} ${WORK_DIR}/bin/test2_timing 3 24 1 3 2000 0 1e-9 1e-14 > ${WORK_DIR}/results/2000ptcls_3peri_80proc_streamlines.txt

# make test_proxy_mats -B -j &&
# mpirun -n 1 --map-by slot:pe=${OMP_NUM_THREADS} ./bin/test2_ptcl_periodize_mats 30 20 > proxy_mat_30_20.txt

# make test_proxy_order -B -j &&
# mpirun -n 1 --map-by slot:pe=${OMP_NUM_THREADS} ./bin/test2_ptcl_periodize_order 4 48 1 25 0 30000 > proxy_order_30_20_debug.txt
# NOTE: above uses regular matrices, works well. 
#       below uses float128 matrices
# mpirun -n 1 --map-by slot:pe=${OMP_NUM_THREADS} ./bin/test2_ptcl_periodize_order 2 16 1 25 0 30000 > proxy_order_2_20_quadreal.txt

# make test1_dispersion -B -j &&
# mpirun -n 6 --map-by numa:pe=${OMP_NUM_THREADS} ./bin/test1_dispersion 400 64 5 > dispersion.txt

