#!/bin/bash

# Set up batch job settings
# #SBATCH --job-name=Convdiv
# #SBATCH --job-name=weak_scaling_3peri_low
#SBATCH --job-name=dash_n_trg
#SBATCH --mail-type=BEGIN,END
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=2
#SBATCH --cpus-per-task=18
#SBATCH --mem-per-cpu=2g
#SBATCH --time=00:15:00
#SBATCH --account=shravan0
#SBATCH --partition=standard
#SBATCH --exclusive

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export KMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export NTASK=$((${SLURM_NNODES}*${SLURM_NTASKS_PER_NODE}))
export NCORES=$((${NTASK}*${OMP_NUM_THREADS}))

mpirun --report-bindings -n 1 --map-by slot:pe=${OMP_NUM_THREADS} ./bin/test2_ptcl_conv 2 32 0 1 100 0 30000 > 100ptcls_1peri_1nodes.txt
mpirun --report-bindings -n 2 -mca coll_hcoll_enable 0 --map-by slot:pe=${OMP_NUM_THREADS} ./bin/test2_ptcl_conv 2 32 0 1 100 0 30000 > 100ptcls_1peri_2nodes.txt
