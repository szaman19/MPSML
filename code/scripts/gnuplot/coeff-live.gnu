#!/usr/local/gnuplot/src/gnuplot

#set terminal qt font ",14" 
set terminal pdfcairo size 8,6
set output "c.pdf"

set multiplot layout 1,2

set border lw 2

set xlabel "Bx"
set format x "%g"
set xtics nomirror

unset key

d0 = "data/output/ferro/4-site-bb-hamil-gs-coeff.dat"
d1 = "data/output/ferro/4-site-pg-hamil-gs-coeff.dat"

set ylabel "C_i^2"

set title "Black Box"
plot for [i=2:17] d0 every 100 u 1:i w p ps 0.5 title "c".i

set title "Physics Guided"
plot for [i=2:17] d1 every 100 u 1:i w p ps 0.5 title "c".i

#pause 2
#reread
