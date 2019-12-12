#!/bin/bash 

#This script generates a gnuplot script on the fly to plot
#single epoch snapshots of Mz. It collates all the PDFS 
#into single "movie" using ghostscript".

module load gnuplot 

usage ()
{
	echo 
	echo "Usage: generate-Mz-movie.sh <chain> <qubits> <sampling>"
	echo 
	exit 1
}

if [ $# -ne 3 ]
then
	usage
fi

PDFS=()

bb=$(find -O3 data/output/$1/$2-qubits -name '*.pkl' | grep $3 | grep bb)
c2=$(find -O3 data/output/$1/$2-qubits -name '*.pkl' | grep $3 | grep c2)
pg=$(find -O3 data/output/$1/$2-qubits -name '*.pkl' | grep $3 | grep pg)
rd=$(find -O3 data/output/$1/$2-qubits -name '*.pkl' | grep $3 | grep rd)

bb=${bb%/opt-bb.pkl}
c2=${c2%/opt-c2.pkl}
pg=${pg%/opt-pg.pkl}
rd=${rd%/opt-rd.pkl}

echo "best bb in $bb"
echo "best c2 in $c2"
echo "best pg in $pg"
echo "best rd in $rd"

rm -f $bb/testing/*.pdf
rm -f $c2/testing/*.pdf
rm -f $pg/testing/*.pdf
rm -f $rd/testing/*.pdf

max_bb=$(ls -1vlAt $bb/testing/*bb* | head -1 | awk -F '=' '{print $2}')
max_c2=$(ls -1vlAt $c2/testing/*c2* | head -1 | awk -F '=' '{print $2}')
max_pg=$(ls -1vlAt $pg/testing/*pg* | head -1 | awk -F '=' '{print $2}')
max_rd=$(ls -1vlAt $rd/testing/*rd* | head -1 | awk -F '=' '{print $2}')

max_bb=${max_bb%.dat}
max_c2=${max_c2%.dat}
max_pg=${max_pg%.dat}
max_rd=${max_rd%.dat}

vals=( $max_bb $max_c2 $max_pg $max_rd )

MAX=${vals[0]}
for v in "${vals[@]}"
do
	((v > max)) && MAX=$v
done

MIN=0
STP=1

if [ "$3" == "single" ]
then 
	for (( epoch = $MIN; epoch <= $MAX; epoch += $STP ))
	do
		cat > tmp.gnu <<- EOF
		set terminal pdfcairo enhanced size 4,3 font ",12"
		set output "data/output/$1/$2-qubits/40k-instances/single-phase/testing/mz-epoch=$epoch.pdf"

		bb = "$bb/testing/avg-bb-epoch=$epoch.dat"
		pg = "$pg/testing/avg-pg-epoch=$epoch.dat"
		c2 = "$c2/testing/avg-c2-epoch=$epoch.dat"
		rd = "$rd/testing/avg-rd-epoch=$epoch.dat"
		mz = "data/output/$1/$2-qubits/40k-instances/single-phase/testing/avg-sz-exact.dat"

		set key outside
		#set style rect fc lt -1 fs solid 0.15 noborder
		set style rect fc rgb "grey" fs pattern 1 bo -1
		set obj rect from graph 0, graph 0 to 0.5, graph 1
		set label "Training" at 0.05,0.2

		set title "Epoch $epoch"
		#set yrange[0:1]

		set xlabel "B_x"
		set ylabel "M_z" rotate by 360

		p mz every 10 u 1:2 w p pt 7 ps 0.5 lc rgb "purple" title "Actual", \\
		  bb every 10 u 1:2 w p pt 7 ps 0.5 lc rgb "dark-grey" title "Black Box",\\
		  pg every 10 u 1:2 w p pt 7 ps 0.5 lc rgb "dark-blue" title "{/Symbol d} Method",\\
		  c2 every 10 u 1:2 w p pt 7 ps 0.5 lc rgb "dark-red" title "C_2 Method", \\
		  rd every 10 u 1:2 w p pt 7 ps 0.5 lc rgb "dark-green" title "Rd Method",
		EOF

		gnuplot tmp.gnu #2> /dev/null
		rm tmp.gnu
		
		PDFS+=( "data/output/$1/$2-qubits/40k-instances/single-phase/testing/mz-epoch=$epoch.pdf" )
	done
else
	for (( epoch = $MIN; epoch <= $MAX; epoch += $STP ))
	do
		cat > tmp.gnu <<- EOF
		set terminal pdfcairo enhanced size 4,3 font ",12"
		set output "data/output/$1/$2-qubits/40k-instances/single-phase/testing/mz-epoch=$epoch.pdf"

		bb = "$bb/testing/avg-bb-epoch=$epoch.dat"
		pg = "$pg/testing/avg-pg-epoch=$epoch.dat"
		c2 = "$c2/testing/avg-c2-epoch=$epoch.dat"
		rd = "$rd/testing/avg-rd-epoch=$epoch.dat"
		mz = "data/output/$1/$2-qubits/40k-instances/single-phase/testing/avg-sz-exact.dat"

		set key outside

		set title "Epoch $epoch"
		#set yrange[0:1]

		set xlabel "B_x"
		set ylabel "M_z" rotate by 360

		p mz every 10 u 1:2 w p pt 7 ps 0.5 lc rgb "purple" title "Actual", \\
		  bb every 10 u 1:2 w p pt 7 ps 0.5 lc rgb "dark-grey" title "Black Box",\\
		  pg every 10 u 1:2 w p pt 7 ps 0.5 lc rgb "dark-blue" title "{/Symbol d} Method",\\
		  c2 every 10 u 1:2 w p pt 7 ps 0.5 lc rgb "dark-red" title "C_2 Method", \\
		  rd every 10 u 1:2 w p pt 7 ps 0.5 lc rgb "dark-green" title "Rd Method",
		EOF

		gnuplot tmp.gnu #2> /dev/null
		rm tmp.gnu
		
		PDFS+=( "data/output/$1/$2-qubits/40k-instances/single-phase/testing/mz-epoch=$epoch.pdf" )
	done
fi

gs -q -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -sOutputFile=$1-$2-site-$3-phase-Mz-movie.pdf ${PDFS[@]}
