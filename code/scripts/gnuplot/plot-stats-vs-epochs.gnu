#!/usr/bin/env gnuplot 


set terminal pdfcairo enhanced size 5,5 font ",11"
set output "learning.pdf"

set multiplot layout 3,2 title "10 Training Points (Randomly sampled, online mode)"

bb = 'bb-2-mean-and-variance-vs-epochs.dat'
c2 = 'c2-2-mean-and-variance-vs-epochs.dat'
rr = 'rr-2-mean-and-variance-vs-epochs.dat'
pv = 'pv-2-mean-and-variance-vs-epochs.dat'
#mx = 'mpv-6-mean-and-variance-vs-epochs.dat'
#sr = 'srr-4-mean-and-variance-vs-epochs.dat'

set logscale x
set xlabel "Epochs"

lw = 2
set border lw 2
#set logscale y
#set yrange[0.999:1]
#set xrange[100:]

set key bottom right
set title "(a) Mean Overlap (2 spins 100 trials)"
p bb u 1:2 w l lw lw lc rgb "red" title "MSE ALL", \
  pv u 1:2 w l lw lw lc rgb "black" title "VAR ALL" , \
  rr u 1:2 w l lw lw lc rgb "purple" title "(MSE+VAR) ALL", \
  c2 u 1:2 w l lw lw lc rgb "blue" title "(MSE+SCH) ALL",\
  #sr u 1:2 w l lw lw lc rgb "dark-green" title "MSE <= 0.5 + sig(VAR) ELSE"

unset key
unset logscale y
set logscale y
#set yrange[0.00000001:0.1]
#set xrange[:]
set title "(b) Overlap Variance (2 spins 100 trials)"
p bb u 1:3 w l lw lw lc rgb "red" title "MSE", \
  pv u 1:3 w l lw lw lc rgb "black" title "VAR", \
  rr u 1:3 w l lw lw lc rgb "purple" title "(MSE+VAR) ALL",\
  c2 u 1:3 w l lw lw lc rgb "blue" title "(MSE+SCH) ALL",\
  #mx u 1:3 w l lw lw lc rgb "blue" title "MSE <= 0.25 + VAR ELSE", \
  #sr u 1:3 w l lw lw lc rgb "dark-green" title "MSE <= 0.5 + sig(VAR) ELSE"

bb = 'bb-4-mean-and-variance-vs-epochs.dat'
c2 = 'c2-4-mean-and-variance-vs-epochs.dat'
rr = 'rr-4-mean-and-variance-vs-epochs.dat'
pv = 'pv-4-mean-and-variance-vs-epochs.dat'

unset logscale y
set key bottom right
set title "(a) Mean Overlap (4 spins 100 trials)"
p bb u 1:2 w l lw lw lc rgb "red" title "MSE ALL", \
  pv u 1:2 w l lw lw lc rgb "black" title "VAR ALL" , \
  rr u 1:2 w l lw lw lc rgb "purple" title "(MSE+VAR) ALL", \
  c2 u 1:2 w l lw lw lc rgb "blue" title "(MSE+SCH) ALL",\
  #sr u 1:2 w l lw lw lc rgb "dark-green" title "MSE <= 0.5 + sig(VAR) ELSE"

unset key
unset logscale y
set logscale y
#set yrange[0.00000001:0.1]
#set xrange[:]
set title "(b) Overlap Variance (4 spins 100 trials)"
p bb u 1:3 w l lw lw lc rgb "red" title "MSE", \
  pv u 1:3 w l lw lw lc rgb "black" title "VAR", \
  rr u 1:3 w l lw lw lc rgb "purple" title "(MSE+VAR) ALL",\
  c2 u 1:3 w l lw lw lc rgb "blue" title "(MSE+SCH) ALL",\
  #mx u 1:3 w l lw lw lc rgb "blue" title "MSE <= 0.25 + VAR ELSE", \
  #sr u 1:3 w l lw lw lc rgb "dark-green" title "MSE <= 0.5 + sig(VAR) ELSE"

bb = 'bb-6-mean-and-variance-vs-epochs.dat'
c2 = 'c2-6-mean-and-variance-vs-epochs.dat'
rr = 'rr-6-mean-and-variance-vs-epochs.dat'
pv = 'pv-6-mean-and-variance-vs-epochs.dat'

unset logscale y
set key bottom right
set title "(a) Mean Overlap (6 spins 100 trials)"
p bb u 1:2 w l lw lw lc rgb "red" title "MSE ALL", \
  pv u 1:2 w l lw lw lc rgb "black" title "VAR ALL" , \
  rr u 1:2 w l lw lw lc rgb "purple" title "(MSE+VAR) ALL", \
  c2 u 1:2 w l lw lw lc rgb "blue" title "(MSE+SCH) ALL",\
  #sr u 1:2 w l lw lw lc rgb "dark-green" title "MSE <= 0.5 + sig(VAR) ELSE"

unset key
unset logscale y
set logscale y
#set yrange[0.00000001:0.1]
#set xrange[:]
set title "(b) Overlap Variance (6 spins 100 trials)"
p bb u 1:3 w l lw lw lc rgb "red" title "MSE", \
  pv u 1:3 w l lw lw lc rgb "black" title "VAR", \
  rr u 1:3 w l lw lw lc rgb "purple" title "(MSE+VAR) ALL",\
  c2 u 1:3 w l lw lw lc rgb "blue" title "(MSE+SCH) ALL",\
  #mx u 1:3 w l lw lw lc rgb "blue" title "MSE <= 0.25 + VAR ELSE", \
  #sr u 1:3 w l lw lw lc rgb "dark-green" title "MSE <= 0.5 + sig(VAR) ELSE"
