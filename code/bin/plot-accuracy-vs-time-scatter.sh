#!/bin/bash 

programname=$0
root=/home/csingh5/Documents/guided-machine-learning/code

function usage {
    echo "usage: $programname [num-qubits]"
}

if [ $# -ne 1 ]
then
	usage
	exit 1
fi

for algo in {bb,c2,rd,pg}
do
	file=$root/data/output/clean-Ising/$1-$algo.txt
	grep "Training time" $file | awk '{print $3}' > tmp1
	grep "Average overlap" $file | awk '{print $3}' > tmp2
	paste tmp1 tmp2 > $root/data/output/clean-Ising/$1-$algo-acc-vs-time.txt
	rm tmp1 tmp2
done


cat > tmp.gnu <<- EOF
	set terminal pdfcairo enhanced size 5,4 font ",10"
	set output "$1-acc-vs-time.pdf"
	set encoding utf8

	set title "Dual phase $1 qubit"

	a1 = '$root/data/output/clean-Ising/$1-bb-acc-vs-time.txt'
	a2 = '$root/data/output/clean-Ising/$1-c2-acc-vs-time.txt'
	a3 = '$root/data/output/clean-Ising/$1-rd-acc-vs-time.txt'
	a4 = '$root/data/output/clean-Ising/$1-pg-acc-vs-time.txt'

	set xlabel "Training Time (s)"
	set ylabel "Average Overlap"
	
	set key outside top horizontal center

	set logscale x

	set yrange[:1]
	#set xrange[0.4:100]

	ps = 0.5

	p a1 w p pt 7 ps ps lc rgb "dark-blue" title "MSE", \
	  a2 w p pt 5 ps ps lc rgb "dark-red" title "MSE + SCH", \
	  a3 w p pt 3 ps ps lc rgb "dark-green" title "MSE + RAND",\
	  a4 w p pt 1 ps ps lc rgb "black" title "PSGD"
EOF

gnuplot tmp.gnu 
rm tmp.gnu 


for algo in {bb,c2,rd,pg}
do
	rm $root/data/output/clean-Ising/$1-$algo-acc-vs-time.txt
done
