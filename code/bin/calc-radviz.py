#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import sys


#number of trials in bin file... 
N = int(sys.argv[1])

with open("wandb-radviz.bin", "rb") as file:
    a = np.fromfile(file)

a = np.reshape(a, (N, -1))

def anchorx(j):
    return np.cos(2*np.pi*j/a.shape[1])

def anchory(j):
    return np.sin(2*np.pi*j/a.shape[1])

max = np.amax(a)
min = np.amin(a)

print(max, min)

# a -= min 
# a /= max - min

for i in range(a.shape[1]):
    a[:,i] = (a[:,i] - np.amin(a[:,i])) / (np.amax(a[:,i]) - np.amin(a[:,i]))

b = np.zeros((N, 2))

for i in range(N):
    denom = np.sum(a[i,:])
    for j in range(a.shape[1]):
        b[i, 0] += a[i,j] * anchorx(j) / denom
        b[i, 1] += a[i,j] * anchory(j) / denom


fig, ax = plt.subplots()

n1 = 2000
x = b[:n1-1,0]; y = b[:n1-1,1]
ax.plot(x, y, label="BB trained")
ax.scatter(x[0], y[0], color="red")
ax.scatter(x[-1], y[-1], color="green")

n2 = 4000
x = b[n1:n2-1,0]; y = b[n1:n2-1,1]
ax.plot(x, y, label="C2 trained")
ax.scatter(x[0], y[0], color="red")
ax.scatter(x[-1], y[-1], color="green")

x = b[n2:,0]; y = b[n2:,1]
ax.plot(x, y, label="Semi trained")
ax.scatter(x[0], y[0], color="red")
ax.scatter(x[-1], y[-1], color="green")

# n2 = N
# ax.plot(b[:n2-1,0], b[:n2-1,1], label="BB trained")
# ax.scatter(b[n1,0], b[n1,1], label="BB trained start", color="red")
# ax.scatter(b[n2-1,0], b[n2-1,1], label="BB trained end", color="green")

# ax.plot(b[1001:2001,0], b[1001:2001,1], label="C2 trained")
# ax.scatter(b[1001,0], b[1001,1], label="C2 trained start", color="red")
# ax.scatter(b[-1,0], b[-1,1], label="C2 trained end", color="green")

# ax.plot(b[2002:,0], b[2002:,1], label="Semi trained")
# ax.scatter(b[2002,0], b[2002,1], label="Semi trained start", color="red")
# ax.scatter(b[-1,0], b[-1,1], label="Semi trained end", color="green")

# ax.scatter(b[300:,0], b[300:,1], label="C2 trained")
# ax.scatter(b[400:,0], b[400:,1], label="Semi-Supervised")
# ax.scatter(b[100:300,0], b[100:300,1], label="untrained")
plt.xlabel("x")
plt.ylabel("y")
plt.title("Radial Visualization in Weight Space over Epochs")
# plt.xlim((-0.08,0.1))
# plt.ylim((-0.08,0.04))

ax.legend()

plt.savefig("tmp.png", dpi=900)
