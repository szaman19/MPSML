#!/usr/bin/env gnuplot 

set terminal pdfcairo enhanced size 3.25,5.5 font ",11"
set output "bars.pdf"

set multiplot layout 4,1

set boxwidth 0.9 relative
set style data histograms
set style fill solid 1.0 border -1

bb1 = "bb-sub-hist.txt"
bb2 = "bb-cri-hist.txt"
bb3 = "bb-sup-hist.txt"

pg1 = "pg-sub-hist.txt"
pg2 = "pg-cri-hist.txt"
pg3 = "pg-sup-hist.txt"

c21 = "c2-sub-hist.txt"
c22 = "c2-cri-hist.txt"
c23 = "c2-sup-hist.txt"

bc = "field-lyapunov-bb.txt"
pc = "field-lyapunov-pg.txt"
cc = "field-lyapunov-c2.txt"

set ylabel "B_x" rotate by 360
set xrange[2:14]
set format y "%2.1f"
set ytics nomirror
set xtics nomirror 
unset key 

set arrow from 2,1 to 14,1 nohead dashtype "-"
set label "Critical Field" at 2.2,1.2
set label "(a)" at 2.5,1.8
set xlabel " "
#set ytics (0,0.5,1.0,1.5,2.0)
set lmargin 8.5
plot bc every 2 u 2:1 w p pt 7 ps 0.05 lc rgb "dark-blue" title "Black Box", \
     pc every 3 u 2:1 w p pt 7 ps 0.05 lc rgb "dark-red" title "{/Symbol d}", \
     cc every 4 u 2:1 w p pt 7 ps 0.05 lc rgb "dark-green" title "C_2"

set lmargin -1
unset label 
unset arrow
set format y "%1.2f"
#set ytics (0,0.4,0.8,0.12)
set key top right 
unset ylabel 
set ylabel "Point Density" 

letter_height = 0.16
actual_height = letter_height - 0.03

set yrange[0:0.18]
set label "(b) Quadratic Cost" at 2.5,letter_height
#set label "Quadratic Cost" at 2.5,actual_height

plot bb1 u 2:($1/12000) with boxes fs transparent solid 0.90 lc rgb "#023858" title "0.0 < B_x < 0.8", \
     bb2 u 2:($1/12000) with boxes fs transparent solid 0.70 lc rgb "#0570b0" title "0.8 < B_x < 1.2", \
     bb3 u 2:($1/12000) with boxes fs transparent solid 0.50 lc rgb "#74a9cf" title "1.2 < B_x < 2.0"

letter_height = 0.1075
actual_height = letter_height - 0.03
arrow_height = 0.05
stability_height = arrow_height - 0.015
set yrange[0:0.12]
unset label 
set arrow from 8.2,arrow_height to 10.2,arrow_height lw 2
set label "Decreasing stability" at 7.5,stability_height
set label "(c) Quadratic + Schrodinger" at 2.5,letter_height
#set label "Quadratic + Schrodinger" at 2.5,actual_height
plot c21 u 2:($1/12000) with boxes fs transparent solid 0.90 lc rgb "#00441b" title "0.0 < B_x < 0.8", \
     c22 u 2:($1/12000) with boxes fs transparent solid 0.70 lc rgb "#238b45" title "0.8 < B_x < 1.2", \
     c23 u 2:($1/12000) with boxes fs transparent solid 0.50 lc rgb "#66c2a4" title "1.2 < B_x < 2.0"

unset arrow 
unset label 
set arrow from 10.2,arrow_height to 8.2,arrow_height lw 2
set xlabel "Lyapunov Estimate"
set label "Increasing stability" at 7.5,stability_height
set label "(d) Perturbation Method" at 2.5,letter_height
#set label "{/Symbol d} Method" at 2.5,actual_height
plot pg1 u 2:($1/12000) with boxes fs transparent solid 0.90 lc rgb "#800026" title "0.0 < B_x < 0.8", \
     pg2 u 2:($1/12000) with boxes fs transparent solid 0.70 lc rgb "#e31a1c" title "0.8 < B_x < 1.2", \
     pg3 u 2:($1/12000) with boxes fs transparent solid 0.50 lc rgb "#fd8d3c" title "1.2 < B_x < 2.0"
