#!/usr/bin/env gnuplot 

set terminal cairolatex standalone pdf size 3.375,3.375 font ",9"
set output "4-site-gs-coeff.tex"

set border lw 2

set multiplot layout 1,2

set border lw 2

set xlabel "Bx"
set xrange[0:1.6]
set xtics 0,0.5,1.5 nomirror

unset key

d0 = "../code/chris/data/output/ferro/4-site-bb-hamil-gs-coeff.dat"
d1 = "../code/chris/data/output/ferro/4-site-pg-hamil-gs-coeff.dat"

set ylabel "Probability"

set label 1 "$|\\uparrow \\uparrow \\uparrow \\uparrow\\rangle$" at 0.5,0.95
set label 2 "$|\\downarrow \\downarrow \\downarrow \\downarrow\\rangle$" at 0.7,0.825



skip = 3
ps = 0.2
lw = 0.5

set title "Black Box"
plot for [i=2:17] d0 every skip u 1:i w p ps ps 



unset label 2
set label 2 "$|\\downarrow \\downarrow \\downarrow \\downarrow\\rangle$" at 0.1,0.2

set title "Physics Guided"
plot for [i=2:17] d1 every skip u 1:i w p ps ps 

