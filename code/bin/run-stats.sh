#!/bin/bash 

#SBATCH --nodelist=europa
#SBATCH --ntasks=1
#SBATCH --output=/dev/null 
#SBATCH --error=/dev/null

programname=$0
root=/home/csingh5/Documents/guided-machine-learning/code

function usage {
    echo "usage: $programname [num-qubits] [loss-type]"
}

if [ $# -ne 2 ]
then
	usage
	exit 1
fi

scontrol update job $SLURM_JOB_ID jobname=$1-$2

sed -i "0,/qubits = xxx/s//qubits = $1/" $root/etc/exhaustive.in
sed -i "0,/loss_type = xxx/s//loss_type = $2/" $root/etc/exhaustive.in

cp $root/etc/exhaustive.in $root/etc/$1-$2.in

sed -i "0,/qubits = $1/s//qubits = xxx/" $root/etc/exhaustive.in
sed -i "0,/loss_type = $2/s//loss_type = xxx/" $root/etc/exhaustive.in

rm -f $root/data/output/clean-Ising/$1-$2.txt

for i in {1..100}
do 
	echo $i
   	phynet $root/etc/$1-$2.in | sed 's/\x1b\[[0-9;]*m//g' >> $root/data/output/clean-Ising/$1-$2.txt
done

rm $root/etc/$1-$2.in

