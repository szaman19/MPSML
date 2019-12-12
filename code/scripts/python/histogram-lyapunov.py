#!/usr/bin/env python

import numpy as np

base_dir = "../../data/output/J-equal-minus-one/single-phase/"
num_epochs = 1000
bx = np.loadtxt(base_dir + "4-site-bb-hamil-gs-lyapunov-epoch=0.dat")[:,0]

def lyapunov(epoch, algorithm):
    return np.loadtxt(base_dir + f"4-site-{algorithm}-hamil-gs-lyapunov-epoch={epoch}.dat")[:,2]

bb = np.empty( len(bx) * num_epochs)
pg = np.empty( len(bx) * num_epochs)
c2 = np.empty( len(bx) * num_epochs)

for epoch in range(0,num_epochs):
    lower = epoch * len(bx)
    upper = epoch * len(bx) + len(bx)
    bb[lower:upper] = lyapunov(epoch, "bb")
    pg[lower:upper] = lyapunov(epoch, "pg")
    c2[lower:upper] = lyapunov(epoch, "c2")
    
bb = bb.reshape( (len(bx), num_epochs), order='F')
pg = pg.reshape( (len(bx), num_epochs), order='F')
c2 = c2.reshape( (len(bx), num_epochs), order='F')

bb = bb.max(axis=1)
pg = pg.max(axis=1)
c2 = c2.max(axis=1)

np.savetxt(base_dir + "4-site-bb-hamil-gs-lyapunov-raw.txt", np.vstack((bx,bb)).T)
np.savetxt(base_dir + "4-site-pg-hamil-gs-lyapunov-raw.txt", np.vstack((bx,pg)).T)
np.savetxt(base_dir + "4-site-c2-hamil-gs-lyapunov-raw.txt", np.vstack((bx,c2)).T)

def calc_hist(data):
    a = np.histogram(data, bins=30)
    b = np.zeros((a[0].size, 2))

    for i in range(a[0].size):
        b[i,0] = a[0][i]
        b[i,1] = a[1][i]

    return b

np.savetxt(base_dir + "4-site-bb-hamil-gs-lyapunov-hist.txt", calc_hist(bb))
np.savetxt(base_dir + "4-site-pg-hamil-gs-lyapunov-hist.txt", calc_hist(pg))
np.savetxt(base_dir + "4-site-c2-hamil-gs-lyapunov-hist.txt", calc_hist(c2))
