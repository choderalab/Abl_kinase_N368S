from simtk.openmm.app import *
from simtk.openmm import *

# set up basic parameters
experiment = 'mt_holo_both_4fs' # setting of the experiment (e.g. different combinations of the CVs)
pdbid = '2GQG' # PDB ID of the system
chain = 'A'
iteration = 250000
work_dir = f'/data/chodera/jiayeguo/projects/cv_selection/sams_simulation/new_trials/{pdbid}_{experiment}_{iteration}'
temperature = 310.15 * unit.kelvin
pressure = 1.0 * unit.atmospheres
ndihedrals = 7 # number of dihedrals we want to restrain
ndistances = 2 # number of distances we want to restrain
targets = [0] # list of dunbrack clusters (sams states) to bias to
coefficient = 1.0 # coefficient for force constant
equi_steps = 10000

prmtop = AmberPrmtopFile('./complex_prep/02.ante.tleap/2GQG.com.wat.leap.prmtop')
inpcrd = AmberInpcrdFile('./complex_prep/02.ante.tleap/2GQG.com.wat.leap.inpcrd')
system = prmtop.createSystem(nonbondedMethod=PME, rigidWater=True, nonbondedCutoff=1*unit.nanometer, constraints = HBonds)
print("Done specifying the system.")
# specify the rest of the context for minimization
integrator = LangevinIntegrator(temperature, 1/unit.picosecond, 0.002*unit.picoseconds)
print("Done specifying integrator.")
platform = Platform.getPlatformByName('CUDA')
print("Done specifying platform.")
platform.setPropertyDefaultValue('Precision', 'mixed')
print("Done setting the precision to mixed.")
minimize = Simulation(prmtop.topology, system, integrator, platform)
print("Done specifying simulation.")
minimize.context.setPositions(inpcrd.positions)
print("Done recording a context for positions.")
if inpcrd.boxVectors is not None:
    minimize.context.setPeriodicBoxVectors(*inpcrd.boxVectors)
minimize.context.setVelocitiesToTemperature(temperature)
print("Done assigning velocities.")

# start minimization
tolerance = 0.1*unit.kilojoules_per_mole/unit.angstroms
print("Done setting tolerance.")
minimize.minimizeEnergy()
print("Done setting energy minimization.")
minimize.reporters.append(StateDataReporter('relax-hydrogens.log', 1000, step=True, temperature=True, potentialEnergy=True, totalEnergy=True, speed=True))
print("Done with energy minimization.")

# start equilibration
#minimize.context.setVelocitiesToTemperature(temperature)
minimize.system.addForce(MonteCarloBarostat(pressure, temperature))
minimize.step(equi_steps)
print("Done 10000 steps of equilibration.")

# output the minimized protein as a shortcut.
positions = minimize.context.getState(getPositions=True).getPositions()
print("Done updating positions.")
PDBFile.writeFile(prmtop.topology,positions,open(f'{pdbid}_chain{chain}_holo_minequi.pdb', 'w'), keepIds=True)
print("Done outputing minimized and equilibrated pdb.")
# clean the context
del minimize.context
