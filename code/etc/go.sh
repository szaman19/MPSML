#!/bin/bash 

outc=/home/csingh5/Documents/guided-machine-learning/code/data/input/clean-Ising
outd=/home/csingh5/Documents/guided-machine-learning/code/data/input/dirty-Ising

rm -vf $outc/*.bin
rm -vf $outd/*.bin

for i in {2,4,6,8,10}
do
	srun --nodelist=europa --job-name=clean-$i datgen.x $i 100000 1 $outc/$i-qubits.bin &
done

for i in {2,4,6,8,10}
do
	srun --nodelist=europa --job-name=dirty-$i datgen.x $i 10000 10 $outd/$i-qubits.bin &
done
