#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import sys


#number of trials in bin file... 
file = sys.argv[1]

a = np.loadtxt(file)

print(a.shape)
