'''
This script estimates the
'''
import sys
import os
import json
import numpy as np
import mdtraj as md
import math
import pymbar

def interact(input_top,input_coord):

    # get topology info from the structure
    topology = md.load(input_top).topology
    table, bonds = topology.to_dataframe()
    atoms = table.values

    # setup toolbar to show progress
    toolbar_width = 50
    sys.stdout.write("[%s]" % (" " * toolbar_width))
    sys.stdout.flush()
    sys.stdout.write("\b" * (toolbar_width+1)) # return to start of line, after '['

    traj = md.load(input_coord,top=input_top)
    hb = {}
    for h in range(9):
        hb[h] = [0] * 501

    for i in range(len(traj)): # for each frame in the input trajectory
        hbonds = md.baker_hubbard(traj.slice(i), periodic=False)
        label = lambda hbond : '%s -- %s' % (traj.topology.atom(hbond[0]), traj.topology.atom(hbond[2]))
        for hbond in hbonds:
            if 'ASN368-ND2 -- ASP363-O' == label(hbond):
                hb[0][i] = 1
            if 'ASN368-ND2 -- ALA380-O' == label(hbond):
                hb[1][i] = 1
            if 'ASN368-ND2 -- ASP381-OD1' == label(hbond):
                hb[2][i] = 1
            if 'ASN368-ND2 -- ASP381-OD2' == label(hbond):
                hb[3][i] = 1
            if 'ASN368-ND2 -- ASP363-OD1' == label(hbond):
                hb[4][i] = 1
            if 'ASN368-ND2 -- ASP363-OD2' == label(hbond):
                hb[5][i] = 1
            if 'ALA365-N -- ASN368-OD1' == label(hbond):
                hb[6][i] = 1
            if 'ARG367-NH1 -- ASN368-OD1' == label(hbond):
                hb[7][i] = 1
            if 'ASN368-ND2 -- HIS361-NE2' == label(hbond):
                hb[8][i] = 1

        # update the bar
        if i % 125 == 0:
            sys.stdout.write("-")
            sys.stdout.flush()
    sys.stdout.write("\n")
    # print results
    #print("hb: ",len(hb[1]))
    # clean up
    del traj
    return hb

###################################################################################
## clean up any existing output file and start fresh
if os.path.isfile("./test_2.txt"):
    os.remove("./test_2.txt")

with open("test_2.txt", "a") as f:
    if os.path.isfile('../trajectories/2.dcd'):
        hb = interact('../trajectories/2.pdb', '../trajectories/2.dcd')
        f.write(json.dumps(hb))
        f.write("\n")
        del hb
    else:
        print(f"Trajectory does not exist. Please double-check.")
    f.close()
