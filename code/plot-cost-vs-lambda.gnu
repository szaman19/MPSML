#!/usr/bin/env gnuplot 

set terminal pdfcairo enhanced size 4,5.5 font ",10"
set output "t.pdf" 

set multiplot layout 4,1

data = "tra.txt"

set xlabel "Lagrange Multiplier"
set ylabel "Quadratic Loss"
set ytics nomirror tc "red"
set y2label "Schrodinger Loss" rotate by 270
set y2tics tc "blue"

#set logscale y
#set logscale y2

unset key

#set xrange[0:1.5]

lw = 2
set title "2 site state 1"
plot data u 1:2 w l lw lw lc rgb "red" title "C_1", \
	 data u 1:6 axes x1y2 w l lw lw lc rgb "blue" title "C_2"

set title "2 site state 2"
plot data u 1:3 w l lw lw lc rgb "red" title "C_1", \
	 data u 1:7 axes x1y2 w l lw lw lc rgb "blue" title "C_2"

set title "2 site state 3"
plot data u 1:4 w l lw lw lc rgb "red" title "C_1", \
	 data u 1:8 axes x1y2 w l lw lw lc rgb "blue" title "C_2"

set title "2 site state 4"
plot data u 1:5 w l lw lw lc rgb "red" title "C_1", \
	 data u 1:9 axes x1y2 w l lw lw lc rgb "blue" title "C_2"

#set title "4 site state 1"
#plot data u 1:2 w l lw lw lc rgb "red" title "C_1", \
     #data u 1:18 axes x1y2 w l lw lw lc rgb "blue" title "C_2"

#set title "4 site state 2"
#plot data u 1:3 w l lw lw lc rgb "red" title "C_1", \
     #data u 1:19 axes x1y2 w l lw lw lc rgb "blue" title "C_2"

#set title "4 site state 3"
#plot data u 1:4 w l lw lw lc rgb "red" title "C_1", \
     #data u 1:20 axes x1y2 w l lw lw lc rgb "blue" title "C_2"

#set title "4 site state 4"
#plot data u 1:5 w l lw lw lc rgb "red" title "C_1", \
     #data u 1:21 axes x1y2 w l lw lw lc rgb "blue" title "C_2"

#set title "4 site state 5"
#plot data u 1:6 w l lw lw lc rgb "red" title "C_1", \
     #data u 1:22 axes x1y2 w l lw lw lc rgb "blue" title "C_2"

#set title "4 site state 6"
#plot data u 1:7 w l lw lw lc rgb "red" title "C_1", \
     #data u 1:23 axes x1y2 w l lw lw lc rgb "blue" title "C_2"

#set title "4 site state 7"
#plot data u 1:8 w l lw lw lc rgb "red" title "C_1", \
     #data u 1:24 axes x1y2 w l lw lw lc rgb "blue" title "C_2"

#set title "4 site state 8"
#plot data u 1:9 w l lw lw lc rgb "red" title "C_1", \
     #data u 1:25 axes x1y2 w l lw lw lc rgb "blue" title "C_2"

#set title "4 site state 9"
#plot data u 1:10 w l lw lw lc rgb "red" title "C_1", \
     #data u 1:23 axes x1y2 w l lw lw lc rgb "blue" title "C_2"

#set title "4 site state 10"
#plot data u 1:11 w l lw lw lc rgb "red" title "C_1", \
     #data u 1:24 axes x1y2 w l lw lw lc rgb "blue" title "C_2"

#set title "4 site state 11"
#plot data u 1:12 w l lw lw lc rgb "red" title "C_1", \
     #data u 1:25 axes x1y2 w l lw lw lc rgb "blue" title "C_2"

#set title "4 site state 12"
#plot data u 1:13 w l lw lw lc rgb "red" title "C_1", \
     #data u 1:26 axes x1y2 w l lw lw lc rgb "blue" title "C_2"

#set title "4 site state 13"
#plot data u 1:14 w l lw lw lc rgb "red" title "C_1", \
     #data u 1:27 axes x1y2 w l lw lw lc rgb "blue" title "C_2"

#set title "4 site state 14"
#plot data u 1:15 w l lw lw lc rgb "red" title "C_1", \
     #data u 1:28 axes x1y2 w l lw lw lc rgb "blue" title "C_2"

#set title "4 site state 15"
#plot data u 1:16 w l lw lw lc rgb "red" title "C_1", \
     #data u 1:29 axes x1y2 w l lw lw lc rgb "blue" title "C_2"

#set title "4 site state 16"
#plot data u 1:17 w l lw lw lc rgb "red" title "C_1", \
     #data u 1:30 axes x1y2 w l lw lw lc rgb "blue" title "C_2"
