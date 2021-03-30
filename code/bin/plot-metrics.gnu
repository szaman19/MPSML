#!/usr/bin/env gnuplot 

set terminal pdfcairo enhanced size 6.25,2 font ",11"
set output "metrics.pdf"

set multiplot layout 1,3

set xlabel "B_x/J"
#set ylabel "B_z"

ps = 0.25

set xrange[0:2]
set xtics 0,0.5,2
#set ytics 0,0.1,0.5

unset key

#set hidden3d
#set dgrid3d 30,30 qnorm 8

#set ztics 0.95,0.01,1
set key bottom right
set title "Overlap"

#sp 'ovr.dat' u 1:2:3 w p pt 7 ps ps lc rgb "blue" title "Predicted"
#sp 'ovr.dat' u 1:2:3 w l lc rgb "blue" title "Predicted"

unset key
set key bottom center
p 1 w p pt 6 ps ps lc rgb "red" title "Ground Truth", \
  'ovr.dat' u 1:3 w p pt 7 ps ps lc rgb "blue" title "Predicted"


#set key top right
set title "Magnetization"

#set ztics 0,0.2,1
   
#sp 'mag.dat' u 1:2:3 w l lc rgb "blue" title "Prediction"

#sp 'mag.dat' u 1:2:4 w p pt 6 ps ps lc rgb "red" title "Ground Truth", \
   #'mag.dat' u 1:2:3 w p pt 7 ps ps lc rgb "blue" title "Prediction"

set key top right
unset key
set arrow from 1,0 to 1,1 nohead dashtype "-"
p 'mag.dat' u 1:3 w p pt 6 ps ps lc rgb "red" title "Ground Truth", \
  'mag.dat' u 1:4 w p pt 7 ps ps lc rgb "blue" title "Predicted"

unset key
#set key bottom right
set title "Entanglement Entropy"
   
#sp 'ent.dat' u 1:2:3 w l lc rgb "blue" title "Prediction", \
   #'ent.dat' u 1:2:4 w l lc rgb "red" title "Exact"
unset arrow 
set arrow from 1,0 to 1,0.8 nohead dashtype "-"
p 'ent.dat' u 1:3  w p pt 6 ps ps lc rgb "red" title "Ground Truth", \
  'ent.dat' u 1:4  w p pt 7 ps ps lc rgb "blue" title "Predicted"
