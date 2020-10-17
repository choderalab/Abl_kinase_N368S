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
## nine 6 x ~500 matrices, each for one H-bond
hb0 = list()
hb1 = list()
hb2 = list()
hb3 = list()
hb4 = list()
hb5 = list()
hb6 = list()
hb7 = list()
hb8 = list()

count = 0 # counting trajectories
hb_count = dict()
keys = list(range(9))
for key in keys:
    hb_count[key] = [0, 0]

N = 0 # calculate the average length of the trajs
trajs = [0,2,4,6,8,10]
for traj in trajs:
    ## read hbond presence from txt files
    readout = ast.literal_eval(open(f'hb_presence_{traj}.txt', 'r').readlines()[0].strip('\n'))
    hb0.append(np.array(readout['0']))
    hb1.append(np.array(readout['1']))
    hb2.append(np.array(readout['2']))
    hb3.append(np.array(readout['3']))
    hb4.append(np.array(readout['4']))
    hb5.append(np.array(readout['5']))
    hb6.append(np.array(readout['6']))
    hb7.append(np.array(readout['7']))
    hb8.append(np.array(readout['8']))
    hb_count[0][0] += np.sum(readout['0'])
    hb_count[1][0] += np.sum(readout['1'])
    hb_count[2][0] += np.sum(readout['2'])
    hb_count[3][0] += np.sum(readout['3'])
    hb_count[4][0] += np.sum(readout['4'])
    hb_count[5][0] += np.sum(readout['5'])
    hb_count[6][0] += np.sum(readout['6'])
    hb_count[7][0] += np.sum(readout['7'])
    hb_count[8][0] += np.sum(readout['8'])
    N += len(readout['0'])
    count += 1
N /= 9
g0 = pymbar.timeseries.statisticalInefficiencyMultiple(hb0)
g1 = pymbar.timeseries.statisticalInefficiencyMultiple(hb1)
g2 = pymbar.timeseries.statisticalInefficiencyMultiple(hb2)
g3 = pymbar.timeseries.statisticalInefficiencyMultiple(hb3)
g4 = pymbar.timeseries.statisticalInefficiencyMultiple(hb4)
g5 = pymbar.timeseries.statisticalInefficiencyMultiple(hb5)
g6 = pymbar.timeseries.statisticalInefficiencyMultiple(hb6)
g7 = pymbar.timeseries.statisticalInefficiencyMultiple(hb7)
g8 = pymbar.timeseries.statisticalInefficiencyMultiple(hb8)
## calculate average presence of the key salt bridges with statistical uncertainty
uncertainty = dict()

# flatten the lists
hb0 = [item for sublist in hb0 for item in sublist]
hb1 = [item for sublist in hb1 for item in sublist]
hb2 = [item for sublist in hb2 for item in sublist]
hb3 = [item for sublist in hb3 for item in sublist]
hb4 = [item for sublist in hb4 for item in sublist]
hb5 = [item for sublist in hb5 for item in sublist]
hb6 = [item for sublist in hb6 for item in sublist]
hb7 = [item for sublist in hb7 for item in sublist]
hb8 = [item for sublist in hb8 for item in sublist]

uncertainty[0] = np.var(hb0) * g0 / N
uncertainty[1] = np.var(hb1) * g1 / N
uncertainty[2] = np.var(hb2) * g2 / N
uncertainty[3] = np.var(hb3) * g3 / N
uncertainty[4] = np.var(hb4) * g4 / N
uncertainty[5] = np.var(hb5) * g5 / N
uncertainty[6] = np.var(hb6) * g6 / N
uncertainty[7] = np.var(hb7) * g7 / N
uncertainty[8] = np.var(hb8) * g8 / N

## final clean up
del hb0, hb1, hb2, hb3, hb4, hb5, hb6, hb7, hb8

## collect and plot the data
import matplotlib
matplotlib.use("TkAgg")
import matplotlib.pyplot as plt

# total number of frames in the combined trajectories 
N_total = 501 + 469 + 501 + 501 + 467 + 501
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
    
ticks = list(range(1, 10))
plt.bar(ticks, y, yerr=error, color='b', alpha=0.7) # well 1
plt.ylim(0.0, 1.0)

plt.show()
