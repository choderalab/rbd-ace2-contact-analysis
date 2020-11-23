from __future__ import division, print_function

import sys
from sys import stdout

import mdtraj as md
import simtk.openmm as mm
import simtk.openmm.app as app
import simtk.unit as unit
from mdtraj.reporters import DCDReporter
from simtk.openmm import XmlSerializer
from simtk.openmm import MonteCarloBarostat
from simtk.openmm.app import CheckpointReporter, PDBFile
from simtk.openmm.app.amberinpcrdfile import AmberInpcrdFile
from simtk.openmm.app.amberprmtopfile import AmberPrmtopFile
from simtk.openmm.app.pdbreporter import PDBReporter
from simtk.openmm.app.statedatareporter import StateDataReporter

# Set parameters
print("Reading parameters...")
pressure = 1.0 * unit.bar
temperature = 310 * unit.kelvin
nonbonded_method = app.PME
constraints = app.HBonds
remove_cm_motion = True
hydrogen_mass = 4.0 * unit.amu  # Using HMR
collision_rate = 1.0 / unit.picoseconds
timestep = 0.004 * unit.picoseconds  # We can use a 4fs timestep with HMR

# Set steps and frequencies
nsteps = 2500000  # 10 ns
report_freq = 100
chk_freq = 500
traj_freq = 1000  # 2500 frames

# Set input file names
amber_prmtop_file = "../02_run_tleap/RBD_ACE2_complex_FullGlycos.prmtop"
amber_inpcrd_file = "../02_run_tleap/RBD_ACE2_complex_FullGlycos.inpcrd"

# Set file names
output_prefix = 'output/'
integrator_xml_filename = "integrator_4fs.xml"
state_xml_filename = "equilibrated_state_10ns.xml"
state_pdb_filename = "equilibrated_state_10ns.pdb"
system_xml_filename = "equilibrated_system_10ns.xml"
checkpoint_filename = "equilibrated_checkpoint_10ns.chk"
traj_output_filename = "equilibrated_traj_10ns.dcd"

# Load the AMBER files
print("Creating OpenMM system from AMBER input files...")
prmtop = AmberPrmtopFile(amber_prmtop_file)
inpcrd = AmberInpcrdFile(amber_inpcrd_file)

system = prmtop.createSystem(
    nonbondedMethod=nonbonded_method,
    constraints=constraints,
    temperature=temperature,
    removeCMMotion=remove_cm_motion,
    hydrogenMass=hydrogen_mass,
)

# Add a barostat to the system
system.addForce(MonteCarloBarostat(pressure, temperature))

# Make and serialize integrator - Langevin dynamics
print("Serializing integrator to %s" % integrator_xml_filename)
integrator = mm.LangevinIntegrator(
    temperature,
    collision_rate,  # Friction coefficient
    timestep
)
with open(output_prefix + integrator_xml_filename, "w") as outfile:
    xml = mm.XmlSerializer.serialize(integrator)
    outfile.write(xml)

# Define the platform to use; CUDA, OpenCL, CPU, or Reference. Or do not specify
# the platform to use the default (fastest) platform
platform = mm.Platform.getPlatformByName("OpenCL")
prop = dict(OpenCLPrecision="mixed")  # Use mixed single/double precision

# Create the Simulation object
sim = app.Simulation(prmtop.topology, system, integrator, platform, prop)

# Set the particle positions
sim.context.setPositions(inpcrd.positions)

# Minimize the energy
print("Minimising energy...")
print(
    "  initial : %8.3f kcal/mol"
    % (
        sim.context.getState(getEnergy=True).getPotentialEnergy()
        / unit.kilocalories_per_mole
    )
)
sim.minimizeEnergy()
print(
    "  final : %8.3f kcal/mol"
    % (
        sim.context.getState(getEnergy=True).getPotentialEnergy()
        / unit.kilocalories_per_mole
    )
)

# set starting velocities:
print("Generating random starting velocities")
sim.context.setVelocitiesToTemperature(temperature * unit.kelvin)

# write limited state information to standard out:
sim.reporters.append(
    StateDataReporter(
        stdout,
        reportInterval=report_freq,
        step=True,
        time=True,
        potentialEnergy=True,
        kineticEnergy=True,
        temperature=True,
        speed=True,
        progress=True,
        remainingTime=True,
        totalSteps=nsteps,
        separator="\t",
    )
)

# Write to checkpoint files regularly:
sim.reporters.append(CheckpointReporter(
    file=output_prefix + checkpoint_filename,
    reportInterval=chk_freq
    )
)

# Write out the trajectory
sim.reporters.append(md.reporters.DCDReporter(
    file=output_prefix + traj_output_filename,
    reportInterval=traj_freq
    )
)

# Run NPT dynamics
print("Running dynamics in the NPT ensemble...")
sim.step(nsteps)

# Save and serialize the final state
print("Serializing state to %s" % state_xml_filename)
state = sim.context.getState(
    getPositions=True,
    getVelocities=True,
    getEnergy=True,
    getForces=True
)
with open(output_prefix + state_xml_filename, "w") as outfile:
    xml = mm.XmlSerializer.serialize(state)
    outfile.write(xml)

# Save the final state as a PDB
print("Saving final state as %s" % state_pdb_filename)
with open(output_prefix + state_pdb_filename, "w") as outfile:
    PDBFile.writeFile(
        sim.topology,
        sim.context.getState(
            getPositions=True,
            enforcePeriodicBox=True).getPositions(),
            file=outfile,
            keepIds=True
    )

# Save and serialize system
print("Serializing system to %s" % system_xml_filename)
system.setDefaultPeriodicBoxVectors(*state.getPeriodicBoxVectors())
with open(output_prefix + system_xml_filename, "w") as outfile:
    xml = mm.XmlSerializer.serialize(system)
    outfile.write(xml)
