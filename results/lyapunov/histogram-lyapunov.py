#!/usr/bin/env python

import numpy as np

# base_dir = "output/J-equal-minus-one/single-phase/"
# num_epochs = 1000
# bx = np.loadtxt(base_dir + "4-site-bb-hamil-gs-lyapunov-epoch=0.dat")[:,0]

# def lyapunov(epoch, algorithm):
    # return np.loadtxt(base_dir + f"4-site-{algorithm}-hamil-gs-lyapunov-epoch={epoch}.dat")[:,2]

# bb = np.empty( len(bx) * num_epochs)
# pg = np.empty( len(bx) * num_epochs)
# c2 = np.empty( len(bx) * num_epochs)

# for epoch in range(0,num_epochs):
    # lower = epoch * len(bx)
    # upper = epoch * len(bx) + len(bx)
    # bb[lower:upper] = lyapunov(epoch, "bb")
    # pg[lower:upper] = lyapunov(epoch, "pg")
    # c2[lower:upper] = lyapunov(epoch, "c2")
    
# bb = bb.reshape( (len(bx), num_epochs), order='F')
# pg = pg.reshape( (len(bx), num_epochs), order='F')
# c2 = c2.reshape( (len(bx), num_epochs), order='F')

# bb = bb.max(axis=1)
# pg = pg.max(axis=1)
# c2 = c2.max(axis=1)

# np.savetxt("bb-cross.txt", np.vstack((bx,bb)).T)
# np.savetxt("pg-cross.txt", np.vstack((bx,pg)).T)
# np.savetxt("c2-cross.txt", np.vstack((bx,c2)).T)

bb = np.loadtxt("field-lyapunov-bb.txt")
pg = np.loadtxt("field-lyapunov-pg.txt")
c2 = np.loadtxt("field-lyapunov-c2.txt")

def calc_hist(data):
    a = np.histogram(data, bins=45)
    b = np.zeros((a[0].size, 2))

    for i in range(a[0].size):
        b[i,0] = a[0][i]
        b[i,1] = a[1][i]

    return b


def run(data, algo):

    c1 = 0;
    c2 = 0;
    c3 = 0;

    i1 = []
    i2 = []
    i3 = []

    for i in range(len(data[:,0])):
        if data[i,0] < 0.8:
            c1 += 1
            i1.append(i)
        elif data[i,0] >= 0.8 and data[i,0] <= 1.2:
            c2 += 1
            i2.append(i)
        else:
            c3 += 1
            i3.append(i)
    
    sub = np.empty( (c1,2) )
    cri = np.empty( (c2,2) )
    sup = np.empty( (c3,2) )

    for i,j in enumerate(i1):
        sub[i,:] = data[j,:]

    for i,j in enumerate(i2):
        cri[i,:] = data[j,:]

    for i,j in enumerate(i3):
        sup[i,:] = data[j,:]

    np.savetxt(algo + "-sub-hist.txt", calc_hist(sub))
    np.savetxt(algo + "-cri-hist.txt", calc_hist(cri))
    np.savetxt(algo + "-sup-hist.txt", calc_hist(sup))


run(bb, "bb")
run(pg, "pg")
run(c2, "c2")
