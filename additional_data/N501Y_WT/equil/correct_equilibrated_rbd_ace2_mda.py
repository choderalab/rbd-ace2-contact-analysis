import MDAnalysis as mda
import argparse

# Read filenames
parser = argparse.ArgumentParser(description='correct equilibrated.pdb')
parser.add_argument('-job_id', type=int, help='job number')
parser.add_argument('-prmtop_file', dest='prmtop_file', type=str, help='the prmtop file created from tleap')
parser.add_argument('-inpcrd_file', dest='inpcrd_file', type=str, help='the inprcd file created from tleap')
parser.add_argument('-ref_file', dest='ref_file', type=str, help='a reference PDB to get the correct unit cell dimensions (e.g. ./output/equilibrated.pdb)')

args = parser.parse_args()

prmtop_file = args.prmtop_file
inpcrd_file = args.inpcrd_file
ref_file = args.ref_file

# Load in the topology from tleap output files 
u = mda.Universe(prmtop_file, inpcrd_file)

#u_dim = mda.Universe(ref_file)
#dimensions = u_dim.dimensions

# RBD
rbd = u.select_atoms("index 0-3007")
new_rbd_resids = [i for i in range(332, 528)]
rbd.residues.resids = new_rbd_resids

rbd_glycans = u.select_atoms("index 3008-3241")
new_rbd_glycan_resids = [i for i in range(1, len(rbd_glycans.residues.resids) + 1)]
rbd_glycans.residues.resids = new_rbd_glycan_resids

# ACE2
ace2 = u.select_atoms("index 3242-14591")
new_ace2_resids = [i for i in range(18, 727)]
ace2.residues.resids = new_ace2_resids

ace2_glycans = u.select_atoms("index 14592-15978")
new_ace2_glycan_resids = [i for i in range(1, len(ace2_glycans.residues.resids) + 1)]
ace2_glycans.residues.resids = new_ace2_glycan_resids

ace2_ions = u.select_atoms("index 15979-15980")
new_ace2_ion_resids = [i for i in range(1, len(ace2_ions.residues.resids) + 1)]
ace2_ions.residues.resids = new_ace2_ion_resids

# Solvent
solvent = u.select_atoms("resname Na+ or resname Cl-")
new_solvent_ion_resids = [i for i in range(1, len(solvent.residues.resids) + 1)]
solvent.residues.resids = new_solvent_ion_resids

# Create the new system by merging each universe
new_system = mda.Merge(rbd, rbd_glycans, ace2, ace2_glycans, ace2_ions, solvent)

# Name each chain
new_system.segments.segids = ['R', 'X', 'C', 'D', 'E', 'Y']

#new_system.dimensions = dimensions

# Write out the new system
new_system.atoms.write(f"./output/{args.job_id}/equilibrated_for_imaging.pdb")
