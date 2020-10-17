'''
Identify H-bonds based on the Baker-Hubbard criterion.

'''

import sys
import os
import numpy as np
import mdtraj as md
import math

def interact(input_top, input_coord):
    '''
    When given a PDB code plus a chain index, this script computes distances as well as features such as intermolecular H-bonding that together define protein-ligand interaction.

    '''
    ## open the output file
    #f = open("hbond.txt", "a")

    ## get topology info from the structure
    topology = md.load(input_top).topology
    table, bonds = topology.to_dataframe()
    atoms = table.values

    ## get the array of atom indices for the calculation of:
    ##       * mean of pairwise distances between each ligand atom and CA of 85 binding pocket residues (an (85*n*2) array where n = # of ligand heavy atoms (usually <= 100) and each row contains indices of the two atoms for each distance)
    #dis = np.zeros(shape=(280, 2), dtype=int, order='C')

    ## calculate the distances for the user-specifed structure (a static structure or an MD trajectory)
    traj = md.load(input_coord,top=input_top)
    print(traj)

    ## setup toolbar to show progress
    toolbar_width = 50
    sys.stdout.write("[%s]" % (" " * toolbar_width))
    sys.stdout.flush()
    sys.stdout.write("\b" * (toolbar_width+1)) # return to start of line, after '['

    # dictionary to save h_bond info
    hb = dict()
    for i in range(len(traj)): # for each frame in the traj
        #traj.slice(i).save_pdb(f'./repr_struct/traj{n}_well1_{i}.pdb') # save structures
        hbonds = md.baker_hubbard(traj.slice(i), periodic=False)
        label = lambda hbond : '%s -- %s' % (traj.topology.atom(hbond[0]), traj.topology.atom(hbond[2]))
        for hbond in hbonds:

            #if 'ASN133-ND2' in label(hbond) or 'ASN133-OD1' in label(hbond): # for WT holo 2HYY structure (renumbered)
            #if 'ASN368-ND2' in label(hbond) or 'ASN368-OD1' in label(hbond): # for WT apo 2HYY/2GQG structure
            #if 'ASN146-ND2' in label(hbond) or 'ASN146-OD1' in label(hbond): # for WT holo 2GQG structure (renumbered)

            #if 'SER133-OG' in label(hbond): # for MT holo 2HYY structure (renumbered)
            if 'SER368-OG' in label(hbond): # for MT apo 2HYY/2GQG structure
            #if 'SER146-OG' in label(hbond): # for MT holo 2GQG structure
                if label(hbond) not in hb.keys():
                    hb[label(hbond)] = 1
                else:
                    hb[label(hbond)] += 1
        del hbonds

        ## update the bar
        if i % 125 == 0:
            sys.stdout.write("-")
            sys.stdout.flush()
    sys.stdout.write("\n")
    print(hb)
    ## clean up
    del traj, hb
    return


## remove existing output and start fresh
#if os.path.isfile('./hbond.txt'):
#    os.remove('./hbond.txt')

if os.path.isfile('./2.dcd'):
    interact('./2.pdb', './2.dcd')
else:
    print("Trajectory does not exist. Please double-check.")
