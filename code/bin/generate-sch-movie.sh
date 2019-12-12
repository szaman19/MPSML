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

for (( epoch = $MIN; epoch <= $MAX; epoch += $STP ))
do
	cat > tmp.gnu <<- EOF
	set terminal pdfcairo enhanced size 4,3 font ",12"
	set output "data/output/$1/$2-qubits/40k-instances/single-phase/testing/mz-epoch=$epoch.pdf"

	bb = "$bb/testing/sch-bb-epoch=$epoch.dat"
	pg = "$pg/testing/sch-pg-epoch=$epoch.dat"
	c2 = "$c2/testing/sch-c2-epoch=$epoch.dat"
	rd = "$rd/testing/sch-rd-epoch=$epoch.dat"

	set key outside
	
	set title "Epoch $epoch"
	set yrange[0:300]

	set xlabel "B_x"
	set ylabel "C_2" rotate by 360

	p bb every 5 u 1:2 w p pt 7 ps 0.25 lc rgb "dark-grey" title "Black Box",\\
	  pg every 5 u 1:2 w p pt 7 ps 0.25 lc rgb "dark-blue" title "{/Symbol d} Method",\\
	  c2 every 5 u 1:2 w p pt 7 ps 0.25 lc rgb "dark-red" title "C_2 Method", \\
	  rd every 5 u 1:2 w p pt 7 ps 0.25 lc rgb "dark-green" title "Rd Method",
	EOF

	gnuplot tmp.gnu 
	rm tmp.gnu
	
	PDFS+=( "data/output/$1/$2-qubits/40k-instances/single-phase/testing/mz-epoch=$epoch.pdf" )
done

gs -q -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -sOutputFile=$1-$2-site-$3-phase-sch-movie.pdf ${PDFS[@]}
