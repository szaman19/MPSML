#!/usr/local/gnuplot/src/gnuplot

set terminal pdfcairo enhanced size 6.75,2.75 font ",12" 
set output "mse-multiplot.pdf"

set multiplot layout 1,3 

set border lw 2

set yrange[0.00001:1.1]
set xlabel "Training Epoch"
set format x "%g"
set logscale x
set format y "10^{%L}"
set logscale y

set key bottom left

bb_mse = "data/output/Ising/4-qubits/40k-instances/single-phase/mse-bb.dat"
pg_mse = "data/output/Ising/4-qubits/40k-instances/single-phase/mse-pg.dat"
c2_mse = "data/output/Ising/4-qubits/40k-instances/single-phase/mse-c2.dat"

#bb_mse = "../../data/output/J-equal-minus-two/single-phase/4-site-gs-mse.dat_bb"
#pg_mse = "../../data/output/J-equal-minus-two/single-phase/4-site-gs-mse.dat_pg"
#c2_mse = "../../data/output/J-equal-minus-two/single-phase/4-site-gs-mse.dat_c2"

set ylabel "MSE (Power of base 10)" 

lw = 4

set style line 1 lw lw lc rgb "dark-red" 
set style line 2 lw lw lc rgb "dark-blue"

set title "(a) Quadratic Cost"
plot bb_mse u 1:2 w l ls 1 title "Training", \
     bb_mse u 1:3 w l ls 2 title "Testing"

set title "(b) Quadratic + Schrodinger Cost"
plot c2_mse u 1:2 w l ls 1 title "Training", \
     c2_mse u 1:3 w l ls 2 title "Testing"

set title "(c) Quadratic Cost + {/Symbol d} modification"
plot pg_mse u 1:2 w l ls 1 title "Training", \
     pg_mse u 1:3 w l ls 2 title "Testing"
