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

# SCRIPT FOR WEAK SCALING (increase problem size and processes together) 
mpirun --report-bindings -n 1 --map-by slot:pe=${OMP_NUM_THREADS} ./bin/test2_ptcl_conv 2 32 0 1 25 0 30000 > 25ptcls_1peri_1nodes.txt
mpirun --report-bindings -n 2 -mca coll_hcoll_enable 0 --map-by slot:pe=${OMP_NUM_THREADS} ./bin/test2_ptcl_conv 2 32 0 1 50 0 30000 > 50ptcls_1peri_2nodes.txt
mpirun --report-bindings -n 4 -mca coll_hcoll_enable 0 --map-by slot:pe=${OMP_NUM_THREADS} ./bin/test2_ptcl_conv 2 32 0 1 100 0 30000 > 100ptcls_1peri_4nodes.txt 
mpirun --report-bindings -n 8 -mca coll_hcoll_enable 0 --map-by slot:pe=${OMP_NUM_THREADS} ./bin/test2_ptcl_conv 2 32 0 1 200 0 30000 > 200ptcls_1peri_8nodes.txt 
mpirun --report-bindings -n 16 -mca coll_hcoll_enable 0 --map-by slot:pe=${OMP_NUM_THREADS} ./bin/test2_ptcl_conv 2 32 0 1 400 0 30000 > 400ptcls_1peri_16nodes.txt 

# SCRIPT FOR WEAK SCALING with 3-periodic, no accuracy check, background pressure flow.
mpirun --report-bindings -n 1 --map-by slot:pe=${OMP_NUM_THREADS} ./bin/test2 2 32 3 25 0 30000 > 25ptcls_3peri_1nodes.txt 
mpirun --report-bindings -n 2 -mca coll_hcoll_enable 0 --map-by slot:pe=${OMP_NUM_THREADS} ./bin/test2 2 32 3 50 0 30000 > 50ptcls_3peri_2nodes.txt 
mpirun --report-bindings -n 4 -mca coll_hcoll_enable 0 --map-by slot:pe=${OMP_NUM_THREADS} ./bin/test2 2 32 3 100 0 30000 > 100ptcls_3peri_4nodes.txt 
mpirun --report-bindings -n 8 -mca coll_hcoll_enable 0 --map-by slot:pe=${OMP_NUM_THREADS} ./bin/test2 2 32 3 200 0 30000 > 200ptcls_3peri_8nodes.txt 