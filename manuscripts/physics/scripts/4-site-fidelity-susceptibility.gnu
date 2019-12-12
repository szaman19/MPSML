#!/usr/bin/env gnuplot 

set terminal cairolatex standalone pdf size 3.375,3.375 font ",9"
set output "4-site-fidelity-susceptibility.tex"

set border lw 2

set xlabel "B$_x$"
set xtics nomirror
set yrange[0:70]
unset ytics
unset key

set title "PRELIMINARY RESULT" 

d0 = "../code/chris/data/output/ferro/4-site-fidelity-susceptibility.txt"

set ylabel "Fidelity Susceptibility"
set label "Critical B$_x$" at 1.02, 35
set label "PGML" at 0.4, 67.5

set label "Black Box" at 0.05,17

set arrow from 0.2,15 to 0.25, 12

#set label "PRELIMINARY" at 0.7,35 rotate by 40

set arrow from 1,0 to 1,70 nohead lw 0.5 dashtype "-"
lw = 4

plot d0 u 1:2 w l lw lw lc rgb "blue" title "Physics Guided training",\
     d0 u 1:3 w l lw lw lc rgb "red" title "Black Box training", \

