#!/bin/bash 

for model in {Ising,XXZ}
do
	for qubits in {2,4,6,8,10}
	do
		for sampling in {single,dual}
		do
			echo "running $model $qubits $sampling"
			generate-Mz-movie.sh $model $qubits $sampling
		done
	done
done
