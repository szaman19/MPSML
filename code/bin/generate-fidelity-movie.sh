#!/bin/bash 

#This script generates a gnuplot script on the fly to plot
#single epoch snapshots of the fidelity. It collates all the PDFS 
#into single "movie" using ghostscript". It needs to be 
#called with the directory containing wavefunction files.

MIN=$2
MAX=$3
STP=$4

PDFS=()

for (( epoch = $MIN; epoch <= $MAX; epoch += $STP ))
do
	cat > tmp.gnu <<- EOF
	set terminal pdfcairo enhanced size 4,3 font ",12"
	set output "../../data/output/$1/magnetization-epoch=$epoch.pdf"

	bb = "../../data/output/$1/4-site-bb-hamil-gs-fidelity-epoch=$epoch.dat"
	pg = "../../data/output/$1/4-site-pg-hamil-gs-fidelity-epoch=$epoch.dat"
	c2 = "../../data/output/$1/4-site-c2-hamil-gs-fidelity-epoch=$epoch.dat"

	set key outside
	#set style rect fc lt -1 fs solid 0.15 noborder
	set style rect fc rgb "grey" fs pattern 1 bo -1
	set obj rect from graph 0, graph 0 to 0.5, graph 1

	set title "Epoch $epoch"
	set yrange[0:5]

	set xlabel "B_x"
	set ylabel "Fidelity"

	p bb every 12 u 1:2 w p pt 7 ps 0.5 lc rgb "dark-grey" title "Black Box",\\
	  pg every 12 u 1:2 w p pt 7 ps 0.5 lc rgb "dark-blue" title "{/Symbol d} Method",\\
	  c2 every 12 u 1:2 w p pt 7 ps 0.5 lc rgb "dark-red" title "C_2 Method", \\
	EOF

	gnuplot tmp.gnu 
	rm tmp.gnu
	
	PDFS+=( "../../data/output/$1/magnetization-epoch=$epoch.pdf" )
done

 #( IFS=$'\n'; echo "${PDFS[*]}" )

gs -q -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -sOutputFile=$5.pdf ${PDFS[@]}

