#!/bin/bash

# Partition             Nodes   S-C-T   Timelimit
# ---------             -----   -----   ---------
# sched_mit_hill        (32)    2-8-1   12:00:00
# sched_mit_raffaele    (32)    2-10-1  12:00:00
# sched_any_quicktest   2       2-8-1   00:15:00
# newnodes              (32)    2-10-1  12:00:00

# Job
#SBATCH --partition=sched_any_quicktest
##SBATCH --partition=sched_engaging_default
##SBATCH --partition=newnodes
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=28
##SBATCH --time=00:15:00
#SBATCH --time=12:00:00
##SBATCH --constraint=centos7
##SBATCH --partition=sched_mit_hill
#SBATCH --partition=sched_mit_pog_r8
#SBATCH -J eps04 # sensible name for the job

# Streams
#SBATCH --output=job_%j.out
#SBATCH --error=job_%j.err

# Issues with parallel filesystem
export HDF5_USE_FILE_LOCKING='FALSE'

# Setup conda and dedalus environment
minicondahome="/home/mkohl/miniconda3"
. $minicondahome/etc/profile.d/conda.sh
##source ~/dedalus/dedalus_modules
conda activate dedalus

# Run scripts
## IF YOU NEED TO RESTART:
#ln -s snapshots/snapshots_s1.h5 restart.h5
mpiexec -n 28 -mca btl_tcp_if_include ib0 python3 main.py
mpiexec -n 28 -mca btl_tcp_if_include ib0 python3 -m dedalus merge_procs snapshots --cleanup
mpiexec -n 28  -mca btl_tcp_if_include ib0 python3 -m dedalus merge_procs time_series --cleanup



#mpiexec python3 plot_slices.py snapshots/*.h5



