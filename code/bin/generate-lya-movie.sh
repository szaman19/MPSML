#!/bin/bash 

#This script generates a gnuplot script on the fly to plot
#single epoch snapshots of Mz. It collates all the PDFS 
#into single "movie" using ghostscript". It needs to be 
#called with the directory containing wavefunction files.

module load gnuplot 

MIN=$4
MAX=$5
STP=$6

PDFS=()

bb=$(find -O3 data/output/$1/$2-qubits -name '*.net' | grep $3 | grep bb)
c2=$(find -O3 data/output/$1/$2-qubits -name '*.net' | grep $3 | grep c2)
pg=$(find -O3 data/output/$1/$2-qubits -name '*.net' | grep $3 | grep pg)
rd=$(find -O3 data/output/$1/$2-qubits -name '*.net' | grep $3 | grep rd)

bb=${bb%/bb-model.net}
c2=${c2%/c2-model.net}
pg=${pg%/pg-model.net}
rd=${rd%/rd-model.net}

write_py_file () 
{
cat > tmp.py << EOF
import numpy as np
bb = np.loadtxt("$bb/testing/lya-bb-epoch=$1.dat")
c2 = np.loadtxt("$c2/testing/lya-c2-epoch=$1.dat")
pg = np.loadtxt("$pg/testing/lya-pg-epoch=$1.dat")
rd = np.loadtxt("$rd/testing/lya-rd-epoch=$1.dat")

def calc_hist(data):
	a = np.histogram(data, bins=45)
	b = np.zeros((a[0].size, 2))
	for i in range(a[0].size):
		b[i,0] = a[0][i]
		b[i,1] = a[1][i]
	return b

def run(data, algo, root):
	c1 = 0; c2 = 0; c3 = 0;
	i1 = []; i2 = []; i3 = []

	for i in range(len(data[:,0])):
		if data[i,0] < 0.8:
			c1 += 1
			i1.append(i)
		elif data[i,0] >= 0.8 and data[i,0] <= 1.2:
			c2 += 1
			i2.append(i)
		else:
			c3 += 1
			i3.append(i)
	
	sub = np.empty( (c1,2) )
	cri = np.empty( (c2,2) )
	sup = np.empty( (c3,2) )

	for i,j in enumerate(i1):
		sub[i,:] = data[j,:]

	for i,j in enumerate(i2):
		cri[i,:] = data[j,:]

	for i,j in enumerate(i3):
		sup[i,:] = data[j,:]

	np.savetxt(root + "lya-sub-hist-" + algo + "-epoch=$1.dat", calc_hist(sub))
	np.savetxt(root + "lya-cri-hist-" + algo + "-epoch=$1.dat", calc_hist(cri))
	np.savetxt(root + "lya-sup-hist-" + algo + "-epoch=$1.dat", calc_hist(sup))

run(bb, "bb", "$bb/testing/")
run(c2, "c2", "$c2/testing/")
run(pg, "pg", "$pg/testing/")
run(rd, "rd", "$rd/testing/")
EOF
}

write_gnu_file ()
{
cat > tmp.gnu << EOF
set terminal pdfcairo enhanced size 3.25,5.5 font ",11"
set output "data/output/$1/$2-qubits/40k-instances/single-phase/testing/lya-epoch=$epoch.pdf"

set multiplot layout 4,1

bb = "$bb/testing/lya-sub-hist-bb-epoch=$3.dat"
pg = "$pg/testing/lya-sub-hist-pg-epoch=$3.dat"
c2 = "$c2/testing/lya-sub-hist-c2-epoch=$3.dat"
rd = "$rd/testing/lya-sub-hist-rd-epoch=$3.dat"

set boxwidth 0.9 relative
set style data histograms
set style fill solid 1.0 border -1

bb1 = "$bb/testing/lya-sub-hist-bb-epoch=$3.dat"
bb2 = "$bb/testing/lya-cri-hist-bb-epoch=$3.dat"
bb3 = "$bb/testing/lya-sup-hist-bb-epoch=$3.dat"

pg1 = "$pg/testing/lya-sub-hist-pg-epoch=$3.dat"
pg2 = "$pg/testing/lya-cri-hist-pg-epoch=$3.dat"
pg3 = "$pg/testing/lya-sup-hist-pg-epoch=$3.dat"

c21 = "$c2/testing/lya-sub-hist-c2-epoch=$3.dat"
c22 = "$c2/testing/lya-cri-hist-c2-epoch=$3.dat"
c23 = "$c2/testing/lya-sup-hist-c2-epoch=$3.dat"

bc = "$bb/testing/lya-bb-epoch=$3.dat"
pc = "$pg/testing/lya-pg-epoch=$3.dat"
cc = "$c2/testing/lya-c2-epoch=$3.dat"
rd = "$rd/testing/lya-rd-epoch=$3.dat"

set ylabel "B_x" rotate by 360
#set xrange[2:14]
#set format y "%2.1f"
set ytics nomirror
set xtics nomirror 
unset key 

#set arrow from 2,1 to 14,1 nohead dashtype "-"
#set label "Critical Field" at 2.2,1.2
#set label "(a)" at 2.5,1.8
set xlabel " "
set lmargin 8.5
plot bc u 2:1 w p pt 7 ps 0.05 lc rgb "dark-blue" title "Black Box", \
	 pc u 2:1 w p pt 7 ps 0.05 lc rgb "dark-red" title "{/Symbol d}", \
	 cc u 2:1 w p pt 7 ps 0.05 lc rgb "dark-green" title "C_2"

set lmargin -1
unset label
unset arrow
set format y "%1.2f"
#set ytics (0,0.4,0.8,0.12)
set key top right
unset ylabel
set ylabel "Point Density"

letter_height = 0.16
actual_height = letter_height - 0.03

#set yrange[0:0.18]
#set label "(b) Quadratic Cost" at 2.5,letter_height
#set label "Quadratic Cost" at 2.5,actual_height

plot bb1 u 2:1 with boxes fs transparent solid 0.90 lc rgb "#023858" title "0.0 < B_x < 0.8", \\
	 bb2 u 2:1 with boxes fs transparent solid 0.70 lc rgb "#0570b0" title "0.8 < B_x < 1.2", \\
	 bb3 u 2:1 with boxes fs transparent solid 0.50 lc rgb "#74a9cf" title "1.2 < B_x < 2.0"

#letter_height = 0.1075
#actual_height = letter_height - 0.03
#arrow_height = 0.05
#stability_height = arrow_height - 0.015
#set yrange[0:0.12]
#unset label
#set arrow from 8.2,arrow_height to 10.2,arrow_height lw 2
#set label "Decreasing stability" at 7.5,stability_height
#set label "(c) Quadratic + Schrodinger" at 2.5,letter_height
#set label "Quadratic + Schrodinger" at 2.5,actual_height
plot c21 u 2:1 with boxes fs transparent solid 0.90 lc rgb "#00441b" title "0.0 < B_x < 0.8", \\
	 c22 u 2:1 with boxes fs transparent solid 0.70 lc rgb "#238b45" title "0.8 < B_x < 1.2", \\
	 c23 u 2:1 with boxes fs transparent solid 0.50 lc rgb "#66c2a4" title "1.2 < B_x < 2.0"

unset arrow
unset label
#set arrow from 10.2,arrow_height to 8.2,arrow_height lw 2
#set xlabel "Lyapunov Estimate"
#set label "Increasing stability" at 7.5,stability_height
#set label "(d) Perturbation Method" at 2.5,letter_height
#set label "{/Symbol d} Method" at 2.5,actual_height
plot pg1 u 2:1 with boxes fs transparent solid 0.90 lc rgb "#800026" title "0.0 < B_x < 0.8", \\
	 pg2 u 2:1 with boxes fs transparent solid 0.70 lc rgb "#e31a1c" title "0.8 < B_x < 1.2", \\
	 pg3 u 2:1 with boxes fs transparent solid 0.50 lc rgb "#fd8d3c" title "1.2 < B_x < 2.0"

EOF
}


for (( epoch = $MIN; epoch <= $MAX; epoch += $STP ))
do
	write_py_file $epoch
	write_gnu_file $1 $2 $epoch

	#gnuplot tmp.gnu 
	#rm tmp.gnu
	
	#PDFS+=( "data/output/$1/$2-qubits/40k-instances/single-phase/testing/lya-epoch=$epoch.pdf" )
done

#gs -q -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -sOutputFile=$1-$2-site-$3-phase-lya-movie.pdf ${PDFS[@]}
