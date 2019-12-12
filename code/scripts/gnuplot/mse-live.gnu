#!/usr/local/gnuplot/src/gnuplot

set terminal qt font ",14" 

#set multiplot layout 2,1 

set border lw 2

set xlabel "Training Epoch"
set format x "%g"
set logscale x
set xtics nomirror
#set key center

d0 = "data/output/single-phase/4-site-gs-mse.dat"

set ylabel "MSE"

#set label "Training dataset" at screen 0.25,0.5
#set label "Testing dataset" at screen 0.25,0.25

#epoch 	 BB train 	 PhyNet train 	 C2 train 	 BB test 	 Phynet test 	 C2 test

plot d0 u 1:2 w l lw 2 lc rgb "red" title "Black Box training", \
	 d0 u 1:4 w l lw 2 lc rgb "blue" title "Physics Guided training",\
	 d0 u 1:3 w l lw 2 lc rgb "orange" title "Black Box test", \
	 d0 u 1:5 w l lw 2 lc rgb "purple" title "Physics Guided test"

pause 1

unset label
reread
