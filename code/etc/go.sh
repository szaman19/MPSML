#!/bin/bash 

orig=/home/csingh5/Documents/guided-machine-learning/code/data/input
dest=/run/media/csingh5/Elements

for dir in {clean-Ising,dirty-Ising}
do
	for i in {2,4,6,8}
	do
		mkdir $dest/$dir
		cp -v $orig/$dir/$i-qubits.bin $dest/$dir
	done 
done

#outc=/home/csingh5/Documents/guided-machine-learning/code/data/input/clean-Ising
#outd=/home/csingh5/Documents/guided-machine-learning/code/data/input/dirty-Ising

#cmd=gendat.x

#rm -vf $outc/*.bin
#rm -vf $outd/*.bin

#for i in {2,4,6,8,10}
#do
	#srun --nodelist=europa --job-name=clean-$i $cmd $i 100000 1 $outc/$i-qubits.bin &
#done

#for i in {2,4,6,8,10}
#do
	#srun --nodelist=europa --job-name=dirty-$i $cmd $i 10000 10 $outd/$i-qubits.bin &
#done
