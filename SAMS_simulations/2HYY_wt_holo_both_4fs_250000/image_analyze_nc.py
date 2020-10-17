from netCDF4 import Dataset
import mdtraj as md
import numpy as np

def print_traj():
    experiment = 'both'
    pdbid = '2HYY'
    iteration = 250000
    # get some info directly from nc file
    ncfile = Dataset('traj_checkpoint.nc', 'r')

    # Get an MDTraj topology from a PDB file that matches the atom ordering in the saved trajectories
    print("Generating trajectory ...")
    traj = md.load(f'{pdbid}_chainA_holo_minequi.pdb')
    mdtraj_topology = traj.topology
    # save unitcell info
    lengths = traj.unitcell_lengths 
    angles = traj.unitcell_angles
    nframes = ncfile.dimensions['iteration'].size

    # select the non water part as reference for imaging the frames
    not_water = mdtraj_topology.select('not water')
    anchor_molecules = [{a for a in mdtraj_topology.atoms if a.index in not_water }]
    # Now get the positions from the trajectory
    replica_index = 0
    positions = ncfile.variables['positions'][:,replica_index,:,:] # gets all frames from the first replica
    # Now create an MDTraj trajectory
    mdtraj_trajectory = md.Trajectory(positions, mdtraj_topology, unitcell_lengths=np.repeat(lengths,nframes,axis=0), unitcell_angles=np.repeat(angles,nframes,axis=0))
    # image the frames based on the anchor molecule
    trajectory_imaged = mdtraj_trajectory.image_molecules(anchor_molecules=anchor_molecules)
    # Write the trajectory in gromacs XTC format
    outfilename = f'{experiment}_{iteration}_traj.dcd'
    trajectory_imaged.save(outfilename)
    return

# print traj
print_traj()
