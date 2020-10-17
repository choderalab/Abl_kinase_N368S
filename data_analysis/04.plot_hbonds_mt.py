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
## five 6 x ~500 matrices, each for one H-bond

hb0 = list()
hb1 = list()
hb2 = list()
hb3 = list()
hb4 = list()

count = 0 # counting trajectories
hb_count = dict()
keys = list(range(5))
for key in keys:
    hb_count[key] = [0, 0]

N = 0 # calculate the average length of the trajs
trajs = [1,3,5,7,9,11]
for traj in trajs:
    ## read hbond presence from txt files
    readout = ast.literal_eval(open(f'hb_presence_{traj}.txt', 'r').readlines()[0].strip('\n'))
    hb0.append(np.array(readout['0']))
    hb1.append(np.array(readout['1']))
    hb2.append(np.array(readout['2']))
    hb3.append(np.array(readout['3']))
    hb4.append(np.array(readout['4']))
    hb_count[0][0] += np.sum(readout['0'])
    hb_count[1][0] += np.sum(readout['1'])
    hb_count[2][0] += np.sum(readout['2'])
    hb_count[3][0] += np.sum(readout['3'])
    hb_count[4][0] += np.sum(readout['4'])
    N += len(readout['0'])
    count += 1

N /= 6
g0 = pymbar.timeseries.statisticalInefficiencyMultiple(hb0)
g1 = pymbar.timeseries.statisticalInefficiencyMultiple(hb1)
g2 = pymbar.timeseries.statisticalInefficiencyMultiple(hb2)
g3 = pymbar.timeseries.statisticalInefficiencyMultiple(hb3)
g4 = pymbar.timeseries.statisticalInefficiencyMultiple(hb4)
## calculate average presence of the key salt bridges with statistical uncertainty
uncertainty = dict()

# flatten the lists
hb0 = [item for sublist in hb0 for item in sublist]
hb1 = [item for sublist in hb1 for item in sublist]
hb2 = [item for sublist in hb2 for item in sublist]
hb3 = [item for sublist in hb3 for item in sublist]
hb4 = [item for sublist in hb4 for item in sublist]

uncertainty[0] = np.var(hb0) * g0 / N
uncertainty[1] = np.var(hb1) * g1 / N
uncertainty[2] = np.var(hb2) * g2 / N
uncertainty[3] = np.var(hb3) * g3 / N
uncertainty[4] = np.var(hb4) * g4 / N

## final clean up
del hb0, hb1, hb2, hb3, hb4

## collect and plot the data
import matplotlib
matplotlib.use("TkAgg")
import matplotlib.pyplot as plt

# total number of frames in the combined trajectories 
N_total = 501 + 464 + 501 + 501 + 501 + 415
## normalize data and populate wtih uncertainty
for key in hb_count.keys():
    hb_count[key] = [float(hb_count[key][0]/N_total), uncertainty[key]]

print("Original order:")
print(hb_count)

## sort the list by hb_presence and find the corresponding uncertainty
y = list()
error = list()
for value in sorted(list(value[0] for value in hb_count.values()), reverse=True):
    y.append(value)
    for key in hb_count.keys():
        if value == hb_count[key][0]:
            error.append(hb_count[key][1])
            break
    
ticks = list(range(1, 6))
plt.bar(ticks, y, yerr=error, color='m', alpha=0.7) # well 1
plt.ylim(0.0, 1.0)

plt.show()
