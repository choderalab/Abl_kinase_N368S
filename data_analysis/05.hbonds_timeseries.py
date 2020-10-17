'''
This script estimates the
'''
import sys
import os
import numpy as np
import mdtraj as md
import math
import pymbar
import ast

## set up the answer lists
## nine 6 x 501 matrices, each for one H-bond
hb2 = list()
hb3 = list()

count = 0 # counting trajectories
trajs = [2,6]
for traj in trajs:
    ## read hbond presence from txt files
    readout = ast.literal_eval(open(f'hb_presence_{traj}.txt', 'r').readlines()[0].strip('\n'))
    hb2.append(np.array(readout['2']))
    hb3.append(np.array(readout['3']))
    count += 1
print(np.sum(hb2[0]), np.sum(hb3[0]))
print(np.sum(hb3[1]), np.sum(hb3[1]))

## collect and plot the data
import matplotlib
matplotlib.use("TkAgg")
import matplotlib.pyplot as plt

f, (ax1, ax2) = plt.subplots(2, 1, sharex=True)
#ax1.scatter(list(range(501)), hb2[1], color='b', alpha=0.4, s=8.0)
#ax2.scatter(list(range(501)), hb3[1], color='r', alpha=0.4, s=8.0)
ax1.scatter(list(range(469)), hb2[0][0:469], color='b', alpha=0.4, s=8.0)
ax2.scatter(list(range(469)), hb3[0][0:469], color='r', alpha=0.4, s=8.0)
ax1.set_yticks([0.0, 1.0])
ax2.set_yticks([0.0, 1.0])
plt.show()

## final clean up
del hb2, hb3
