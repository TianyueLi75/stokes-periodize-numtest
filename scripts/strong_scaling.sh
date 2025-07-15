#!/bin/bash

# Set up batch job settings
# #SBATCH --job-name=Convdiv
# #SBATCH --job-name=weak_scaling_3peri_low
#SBATCH --job-name=dash_n_trg
#SBATCH --mail-type=BEGIN,END
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=18
#SBATCH --mem-per-cpu=5g
#SBATCH --time=00:15:00
#SBATCH --account=shravan0
#SBATCH --partition=standard
#SBATCH --exclusive

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export KMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export NTASK=$((${SLURM_NNODES}*${SLURM_NTASKS_PER_NODE}))
export NCORES=$((${NTASK}*${OMP_NUM_THREADS}))



# SCRIPT FOR 500 PTCLS STRONG SCALING (fix problem, increase processes) -- TODO: maybe 3-periodic?
mpirun --report-bindings -n ${NTASK} -mca coll_hcoll_enable 0 --map-by slot:pe=${OMP_NUM_THREADS} ./bin/test2_ptcl_conv 8 64 0 1 500 0 30000 > 500ptcls_1peri_8_64_2proc.txt 
mpirun --report-bindings -n ${NTASK} -mca coll_hcoll_enable 0 --map-by slot:pe=${OMP_NUM_THREADS} ./bin/test2_ptcl_conv 8 64 0 1 500 0 30000 > 500ptcls_1peri_8_64_4proc.txt 
mpirun --report-bindings -n ${NTASK} -mca coll_hcoll_enable 0 --map-by slot:pe=${OMP_NUM_THREADS} ./bin/test2_ptcl_conv 8 64 0 1 500 0 30000 > 500ptcls_1peri_8_64_8proc.txt 
mpirun --report-bindings -n ${NTASK} -mca coll_hcoll_enable 0 --map-by slot:pe=${OMP_NUM_THREADS} ./bin/test2_ptcl_conv 8 64 0 1 500 0 30000 > 500ptcls_1peri_8_64_16proc.txt 
