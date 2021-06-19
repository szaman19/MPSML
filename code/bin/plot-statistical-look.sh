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
	grep "Epochs used" $file | awk '{print $3}' > tmp3
	paste tmp1 tmp2 > $root/data/output/clean-Ising/$1-$algo-acc-vs-time.txt
	paste tmp1 tmp3 > $root/data/output/clean-Ising/$1-$algo-epo-vs-time.txt
	rm tmp1 tmp2 tmp3
done


cat > tmp.gnu <<- EOF
set terminal pdfcairo enhanced size 6.25,3 font ",10"
set output "$1-stats-plot.pdf"
set encoding utf8

set multiplot layout 1,2 title "Dual phase $1 qubit"

a1 = '$root/data/output/clean-Ising/$1-bb-epo-vs-time.txt'
a2 = '$root/data/output/clean-Ising/$1-c2-epo-vs-time.txt'
a3 = '$root/data/output/clean-Ising/$1-rd-epo-vs-time.txt'
a4 = '$root/data/output/clean-Ising/$1-pg-epo-vs-time.txt'

b1 = '$root/data/output/clean-Ising/$1-bb-acc-vs-time.txt'
b2 = '$root/data/output/clean-Ising/$1-c2-acc-vs-time.txt'
b3 = '$root/data/output/clean-Ising/$1-rd-acc-vs-time.txt'
b4 = '$root/data/output/clean-Ising/$1-pg-acc-vs-time.txt'

set xlabel "Training Time (s)"
set ylabel "Epochs used (thousands)"

#set key outside top horizontal center
set key opaque

#set format y "10^{%T}"
#set format y "%.0s{/Symbol \267}10^{%T}"
set format y "%.0s"
#set logscale x

set yrange[0:9999]
#set xrange[0:2000]

ps = 0.5
lw = 2

FIT_LIMIT = 1e-8
f1(x) = m1*x + y0
f2(x) = m2*x + y0
f3(x) = m3*x + y0
f4(x) = m4*x + y0

y0 = 0
fit f1(x) a1 u 1:2 via m1
fit f2(x) a2 u 1:2 via m2
fit f3(x) a3 u 1:2 via m3
fit f4(x) a4 u 1:2 via m4

p a1 w p pt 7 ps ps lc rgb "dark-blue" notitle, \
  f1(x) w l lw lw lc rgb "dark-blue" title "MSE", \
  a2 w p pt 5 ps ps lc rgb "dark-red" notitle, \
  f2(x) w l lw lw lc rgb "dark-red" title "MSE + SCH", \
  a3 w p pt 3 ps ps lc rgb "dark-green" notitle, \
  f3(x) w l lw lw lc rgb "dark-green" title "MSE + RAND", \
  a4 w p pt 1 ps ps lc rgb "black" notitle, \
  f4(x) w l lw lw lc rgb "black" title "PSGD"

######################################################

set xlabel "Training Time (s)"
set ylabel "Average Overlap"

#set key outside top horizontal center

set logscale x
set format y 

set yrange[0.98:1]
#set xrange[0.6:600]

ps = 0.5

p b1 w p pt 7 ps ps lc rgb "dark-blue" title "MSE", \
  b2 w p pt 5 ps ps lc rgb "dark-red" title "MSE + SCH", \
  b3 w p pt 3 ps ps lc rgb "dark-green" title "MSE + RAND",\
  b4 w p pt 1 ps ps lc rgb "black" title "PSGD"
EOF

gnuplot tmp.gnu #2> /dev/null
rm tmp.gnu 
rm fit.log


for algo in {bb,c2,rd,pg}
do
	rm $root/data/output/clean-Ising/$1-$algo-epo-vs-time.txt
	rm $root/data/output/clean-Ising/$1-$algo-acc-vs-time.txt
done
