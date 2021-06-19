#!/bin/bash 

#SBATCH --nodelist=europa

MIN=0.0
MAX=2.0
STP=0.001

rm -vf lam.txt

for lambda in $(seq $MIN $STP $MAX)
do
	sed -i "s/lagrange_multiplier = xx/lagrange_multiplier = $lambda/g" etc/phynet.in
	phynet | tail -1 >> lam.txt
	sed -i "s/lagrange_multiplier = $lambda/lagrange_multiplier = xx/g" etc/phynet.in
	tail -1  lam.txt
done
