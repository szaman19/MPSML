#!/usr/local/gnuplot/src/gnuplot

#set terminal pngcairo size 1000,600 font ",9"
#set output "phase-diagram.png"

d0 = "data/output/whole_diagram_revised/4-site-bb-hamil-un-coeff.dat"
d1 = "data/output/whole_diagram_revised/4-site-pg-hamil-un-coeff.dat"

set xlabel "B_x"
set ylabel "B_z"

set zrange[0:1]
set cbrange[0:1]

unset key 
set title "Predicted phase diagram"

set hidden3d 
set dgrid 50,50 qnorm 2 
set pm3d at b

set view 65,315

sp d1 u 1:2:($3*$3) w l

pause 1
reread
