#!/usr/bin/env gnuplot 

set terminal cairolatex standalone pdf size 3.375,3.375 font ",9"
set output "4-site-gs-mse.tex"

set border 3 lw 2
unset key

set xlabel "Training Epoch"
set logscale x
set xtics nomirror
set ytics nomirror
set xrange[1:600]
set yrange[0:1.3]

d0 = "../code/chris/data/output/ferro/4-site-gs-mse.dat"

set ylabel "MSE"

set label "Physics guided training" at 2,0.75
set arrow from 5,0.725 to 2.1, 0.675

set label "Physics guided test" at 4,0.2
set arrow from 10,0.175 to 10,0.1

set label "Black box training" at 2,1.25
set arrow from 6,1.21 to 6,1.11

set label "Black box test" at 90,0.35
set arrow from 250,0.3 to 250, 0.1

lw = 4

plot d0 u 1:2 w l lw lw lc rgb "red" title "Black Box training", \
	 d0 u 1:3 w l lw lw lc rgb "blue" title "Physics Guided training",\
	 d0 u 1:4 w l lw lw lc rgb "orange" title "Black Box test", \
	 d0 u 1:5 w l lw lw lc rgb "purple" title "Physics Guided test"

