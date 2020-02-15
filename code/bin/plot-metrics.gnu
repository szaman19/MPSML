#!/usr/bin/env gnuplot 

set terminal pdfcairo enhanced size 6.25,2.5 font ",11"
set output "metrics.pdf"

set multiplot layout 1,3

set xlabel "Transverse Field (B_x)"

ps = 0.3


set xrange[0:2]
set xtics 0,0.5,2
set title "Overlap"
set key bottom right
p 1 w p pt 6 ps ps lc rgb "red" title "Ground Truth", \
  'ovr.dat' w p pt 7 ps ps lc rgb "blue" title "Predicted"

set key top right
set title "Magnetization"
p 'mag.dat' u 1:2 w p pt 6 ps ps lc rgb "red" title "Ground Truth", \
  'mag.dat' u 1:3 w p pt 7 ps ps lc rgb "blue" title "Prediction"

set key bottom right
set title "Entanglement Entropy"
p 'ent.dat' u 1:2 w p pt 6 ps ps lc rgb "red" title "Ground Truth", \
  'ent.dat' u 1:3 w p pt 7 ps ps lc rgb "blue" title "Prediction"
