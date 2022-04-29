#!/bin/bash
#SBATCH --job-name=Hamiltonian_Solver_agrace
#SBATCH --mail-type=END        # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=agrace2@binghamton.edu
#SBATCH --cpus-per-task=1            # Number of cores per MPI task 
#SBATCH -N 1            # Maximum number of nodes to be allocated
#SBATCH -p RM-shared
#SBATCH --ntasks-per-node=3
#SBATCH --acctg-freq=task=1
#SBATCH --distribution=cyclic:cyclic # Distribute tasks cyclically first among nodes and then among sockets within a node
#SBATCH --mem=6000MB          # Memory (i.e. RAM) per processor
#SBATCH --time=00:05:00              # Wall time limit (days-hrs:min:sec)
#SBATCH --output=slurm_logs/mpi_test_%j.log     # Path to the standard output and error files relative to the working directory

echo "Date              = $(date)"
echo "Hostname          = $(hostname -s)"
echo "Working Directory = $(pwd)"
echo ""
echo "Number of Nodes Allocated      = $SLURM_JOB_NUM_NODES"
echo "Number of Tasks Allocated      = $SLURM_NTASKS"
echo "Number of Cores/Task Allocated = $SLURM_CPUS_PER_TASK"

#Change arguments here
srun --mpi=pmi2 ./matgen 16 range 1 5 0.0 5.0 5 0.0 5.0 --run-performance-metrics --use-petsc-only-methods