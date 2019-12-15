#!/bin/bash 

rm -vf data/input/clean-Ising/*.bin
rm -vf data/input/dirty-Ising/*.bin

srun --nodelist=europa time datgen.x 2 100000 1 data/input/clean-Ising/2-qubits.bin &
srun --nodelist=europa time datgen.x 4 100000 1 data/input/clean-Ising/4-qubits.bin &
srun --nodelist=europa time datgen.x 6 100000 1 data/input/clean-Ising/6-qubits.bin &
srun --nodelist=europa time datgen.x 8 100000 1 data/input/clean-Ising/8-qubits.bin &
srun --nodelist=europa time datgen.x 10 100000 1 data/input/clean-Ising/10-qubits.bin &

srun --nodelist=europa time datgen.x 2 10000 10 data/input/dirty-Ising/2-qubits.bin &
srun --nodelist=europa time datgen.x 4 10000 10 data/input/dirty-Ising/4-qubits.bin &
srun --nodelist=europa time datgen.x 6 10000 10 data/input/dirty-Ising/6-qubits.bin &
srun --nodelist=europa time datgen.x 8 10000 10 data/input/dirty-Ising/8-qubits.bin &
srun --nodelist=europa time datgen.x 10 10000 10 data/input/dirty-Ising/10-qubits.bin &
