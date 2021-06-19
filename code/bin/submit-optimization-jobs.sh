#!/bin/bash 

rm -rf data/output
xargs mkdir -p < etc/dirs.txt

for model in {Ising,}
do
	for qubits in {2,4,6,8}
	do
		for sampling in {single,dual}
		do
			for algorithm in {bb,pg,c2,rd}
			do
				cat > tmp.sh <<- EOF
				#!/bin/bash 

				#SBATCH --ntasks=1
				#SBATCH --exclude=ganymede
				#SBATCH --nodes=1
				#SBATCH --job-name=$model-$qubits-$sampling-$algorithm
				##SBATCH --mail-type=END
				##SBATCH --mail-user=csingh5@binghamton.edu
				#SBATCH --output=data/output/slurm-files/$model-$qubits-$sampling-$algorithm-slurm-%j.out

				module load psxe-2019
				source activate qm
				export KMP_INIT_AT_FORK=FALSE

				gaussian-bayes.py $model $qubits $sampling $algorithm 
				EOF

				sbatch tmp.sh 
				rm tmp.sh 
			done 
		done 
	done 
done
