#!/bin/bash

# Set up batch job settings
#SBATCH --job-name=jobname
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=2
#SBATCH --cpus-per-task=32
#SBATCH --time=03:35:00
#SBATCH --partition=gen

WORK_DIR=~/stokes-periodize-numtest
source ${WORK_DIR}/pvfmm_modules

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export KMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export NTASK=$((${SLURM_NNODES}*${SLURM_NTASKS_PER_NODE}))
export NCORES=$((${NTASK}*${OMP_NUM_THREADS}))

cd ${WORK_DIR}

# Test periodicity (sanity check)
make test_peri -j &&
mpirun -n 1 --map-by numa:pe=${OMP_NUM_THREADS} ${WORK_DIR}/bin/test_periodicity 3 24 1 1 1e-7 1e-8 > ${WORK_DIR}/out/1ptcls_1peri_periodicity.txt
mpirun -n 1 --map-by numa:pe=${OMP_NUM_THREADS} ${WORK_DIR}/bin/test_periodicity 3 24 2 3 1e-7 1e-8 > ${WORK_DIR}/out/3ptcls_2peri_periodicity.txt
mpirun -n 1 --map-by numa:pe=${OMP_NUM_THREADS} ${WORK_DIR}/bin/test_periodicity 3 24 3 1 1e-7 1e-8 > ${WORK_DIR}/out/1ptcls_3peri_periodicity.txt

# Test self-convergence (Figure 6)
make test_selfconv -j &&
# Singly periodic cylinder-bound flow past a sphere
mpirun -n 1 --map-by numa:pe=${OMP_NUM_THREADS} ${WORK_DIR}/bin/test_convergence 0 1 1 0 > ${WORK_DIR}/test_convergence_channel_1sphere.txt 
# Doubly periodic plane-bound flow past 3 spheres
mpirun -n 1 --map-by numa:pe=${OMP_NUM_THREADS} ${WORK_DIR}/bin/test_convergence 1 2 3 0 > ${WORK_DIR}/test_convergence_plane_3spheres.txt
# Triply periodic 25-sphere-suspension
mpirun -n 1 --map-by numa:pe=${OMP_NUM_THREADS} ${WORK_DIR}/bin/test_convergence 2 3 25 0 > ${WORK_DIR}/test_convergence_25spheres_3peri.txt

# Manufactured solutions test (Figure 7 and Table 3)
make test_manufactured_soln -j &&
for n in {1..8}; do
    for m in {4..80..4}; do
        echo "n = $n, m=$m"
        mpirun -n 1 --map-by numa:pe=${OMP_NUM_THREADS} ./bin/test_manufactured_soln $n $m 25 0 60000 >> "out/Nelem_Nf_grid.txt"
    done
done

# Streamlines (Figure 9)
make timing -j &&
mpirun -n 1 --map-by numa:pe=${OMP_NUM_THREADS} ${WORK_DIR}/bin/timing 4 32 1 3 25 0 1e-9 1e-14 > ${WORK_DIR}/25ptcls_3peri_1proc_streamlines.txt
mpirun -n 16 --map-by numa:pe=${OMP_NUM_THREADS} ${WORK_DIR}/bin/timing 4 32 1 3 400 0 1e-9 1e-14 > ${WORK_DIR}/400ptcls_3peri_16proc_streamlines.txt
mpirun -n 80 --map-by numa:pe=${OMP_NUM_THREADS} ${WORK_DIR}/bin/timing 4 32 1 3 2000 0 1e-9 1e-14 > ${WORK_DIR}/2000ptcls_3peri_80proc_streamlines.txt

# Precompute time (Figure 8)
make timing_precomp -j && 
runs=(
    "1e-2:2"
    "1e-4:4"
    "1e-6:6"
    "1e-8:8"
    "1e-9:10"
    "1e-11:12"
    "1e-13:14"
    "1e-15:16"
    "1e-17:18"
)
# Loop through each pair
for run in "${runs[@]}"; do
    # Extract the value before the colon
    val="${run%%:*}"
    # Extract the suffix number after the colon
    suffix="${run##*:}"

    # Execute the mpirun command
    mpirun -n 1 --report-bindings --map-by numa:pe="${OMP_NUM_THREADS}" \
        "${WORK_DIR}/bin/precompute_time" "${val}" 1.0 0.0 \
        > "${WORK_DIR}/Precomp_timing_FxU_m${suffix}.txt"
done

# Periodization overhead (Table 4)
make timing_periodization -j && 
for np in 25 50 100 200; do
    mpirun -n 1 --report-bindings --map-by numa:pe=${OMP_NUM_THREADS} ${WORK_DIR}/bin/periodization_time 6 32 $np 1e-8 1e-11 > ${WORK_DIR}/periodization_time_Np6Nf32N${np}.txt
done

# Visualization examples (Figure 1a,b)
make examples -j &&
mpirun -n ${NTASK} --map-by numa:pe=${OMP_NUM_THREADS} ${WORK_DIR}/bin/examples 0 > ${WORK_DIR}/channel_example.txt
mpirun -n ${NTASK} --map-by numa:pe=${OMP_NUM_THREADS} ${WORK_DIR}/bin/examples 1 > ${WORK_DIR}/plane_example.txt
