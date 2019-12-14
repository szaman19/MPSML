#!/bin/bash 

srun --nodelist=europa time datgen.x 4 100000 1 data/input/clean-Ising/4-qubits.bin &
srun --nodelist=europa time datgen.x 6 100000 1 data/input/clean-Ising/6-qubits.bin &
srun --nodelist=europa time datgen.x 8 100000 1 data/input/clean-Ising/8-qubits.bin &

srun --nodelist=europa time datgen.x 4 10000 10 data/input/dirty-Ising/4-qubits.bin &
srun --nodelist=europa time datgen.x 6 10000 10 data/input/dirty-Ising/6-qubits.bin &
srun --nodelist=europa time datgen.x 8 10000 10 data/input/dirty-Ising/8-qubits.bin &
