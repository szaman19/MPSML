#!/bin/bash
#SBATCH --job-name=PSOHG
#SBATCH --mail-type=END,FAIL         # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=agrace2@binghamton.edu
#SBATCH --ntasks=12                   # Number of MPI tasks (i.e. processes)
#SBATCH --cpus-per-task=1            # Number of cores per MPI task 
#SBATCH --nodes=4                    # Maximum number of nodes to be allocated
#SBATCH --ntasks-per-node=6          # Maximum number of tasks on each node
#SBATCH --ntasks-per-socket=6        # Maximum number of tasks on each socket
#SBATCH --distribution=cyclic:cyclic # Distribute tasks cyclically first among nodes and then among sockets within a node
#SBATCH --mem-per-cpu=4096mb          # Memory (i.e. RAM) per processor
#SBATCH --time=1-12:00:00              # Wall time limit (days-hrs:min:sec)
#SBATCH --output=slurm_logs/mpi_test_%j.log     # Path to the standard output and error files relative to the working directory

echo "Date              = $(date)"
echo "Hostname          = $(hostname -s)"
echo "Working Directory = $(pwd)"
echo ""
echo "Number of Nodes Allocated      = $SLURM_JOB_NUM_NODES"
echo "Number of Tasks Allocated      = $SLURM_NTASKS"
echo "Number of Cores/Task Allocated = $SLURM_CPUS_PER_TASK"

srun --mpi=pmix_v2 psohg 2 1 1 1 --verbose