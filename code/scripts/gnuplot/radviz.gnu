#!/usr/bin/env gnuplot 

set terminal pdfcairo enhanced size 4,3
set output "radviz.pdf"

set xrange[-1.2:1.2]
set yrange[-1.2:1.2]

#set key center

p 'anchors.dat' w p pt 9 ps 1 lc rgb "dark-blue" title "Anchors", \
  'rad.dat' u 1:2 w p pt 4 ps 0.25 lc rgb "dark-green" title "Truth", \
  'rad.dat' u 3:4 w p pt 5 ps 0.25 lc rgb "dark-red" title "Net"
