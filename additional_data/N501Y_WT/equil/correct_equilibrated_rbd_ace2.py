from simtk.openmm import app, unit
import numpy as np
import argparse
import subprocess

# Read filename
parser = argparse.ArgumentParser(description='correct equilibrated.pdb')
parser.add_argument('name', type=str, help='name of file for which to correct')
args = parser.parse_args()
pdb = app.PDBFile(args.name)

# Get old topology and positions
old_topology = pdb.getTopology()
old_positions = pdb.positions

# Create new topology and positions
new_topology = app.Topology()
new_topology.setPeriodicBoxVectors(old_topology.getPeriodicBoxVectors())
positions_R, positions_X, positions_C, positions_D, positions_E, positions_Y = list(), list(), list(), list(), list(), list()

# Create new chains
new_chain_R = new_topology.addChain(id="R")
new_chain_X = new_topology.addChain(id="X")
new_chain_C = new_topology.addChain(id="C")
new_chain_D = new_topology.addChain(id="D")
new_chain_E = new_topology.addChain(id="E")
new_chain_Y = new_topology.addChain(id="Y")

# Specify the starting residue ids for each chain
d_current_start = {"C": 18, "E": 1, "R": 332, "X": 1, "D": 1}

# Copy residues and atoms to new topology and create split into multiple chains. 
# Also rename residues based on d_current_start 
d_old_to_new = {} # Key: atom in old topology, Value: atom in new topology 
for res in old_topology.residues():  
    residue_id = int(res.id)
    if res.name not in ['HOH', 'Na+', 'Cl-']:
        if residue_id <= 196:
            new_res = new_topology.addResidue(res.name, new_chain_R, id=str(d_current_start["R"]), insertionCode=res.insertionCode)
            for atom in res.atoms():
                new_atom = new_topology.addAtom(atom.name, atom.element, new_res)
                d_old_to_new[atom] = new_atom
                positions_R.append(old_positions[atom.index])
            d_current_start["R"] += 1
        elif residue_id >= 197 and residue_id <= 206:
            new_res = new_topology.addResidue(res.name, new_chain_X, id=str(d_current_start["X"]), insertionCode=res.insertionCode)
            for atom in res.atoms():
                new_atom = new_topology.addAtom(atom.name, atom.element, new_res)
                d_old_to_new[atom] = new_atom
                positions_X.append(old_positions[atom.index])
            d_current_start["X"] += 1
        elif residue_id >= 207 and residue_id <= 915:
            new_res = new_topology.addResidue(res.name, new_chain_C, id=str(d_current_start["C"]), insertionCode=res.insertionCode)
            for atom in res.atoms():
                new_atom = new_topology.addAtom(atom.name, atom.element, new_res)
                d_old_to_new[atom] = new_atom
                positions_C.append(old_positions[atom.index])
            d_current_start["C"] += 1
        elif residue_id >= 916 and residue_id <= 973:
            new_res = new_topology.addResidue(res.name, new_chain_D, id=str(d_current_start["D"]), insertionCode=res.insertionCode)
            for atom in res.atoms():
                new_atom = new_topology.addAtom(atom.name, atom.element, new_res)
                d_old_to_new[atom] = new_atom
                positions_D.append(old_positions[atom.index])
            d_current_start["D"] += 1
        elif residue_id >= 974 and residue_id <= 975:
            new_res = new_topology.addResidue(res.name, new_chain_E, id=str(d_current_start["E"]), insertionCode=res.insertionCode)
            for atom in res.atoms():
                new_atom = new_topology.addAtom(atom.name, atom.element, new_res)
                d_old_to_new[atom] = new_atom
                positions_E.append(old_positions[atom.index])
            d_current_start["E"] += 1
    else:
        new_res = new_topology.addResidue(res.name, new_chain_Y, id=res.id, insertionCode=res.insertionCode)
        for atom in res.atoms():
            new_atom = new_topology.addAtom(atom.name, atom.element, new_res)
            d_old_to_new[atom] = new_atom
            positions_Y.append(old_positions[atom.index])
    
# Copy bonds to new topology
for bond in old_topology.bonds():
    atom_1 = bond[0]
    atom_2 = bond[1]
    atom_1_new = d_old_to_new[atom_1]
    atom_2_new = d_old_to_new[atom_2]
    new_topology.addBond(atom_1_new, atom_2_new)

# Combine positions into list with the right order of the new positions
new_positions = positions_R + positions_X + positions_C + positions_D + positions_E + positions_Y

# Correct new positions to use one Quantity object
new_positions_corrected = unit.quantity.Quantity(value = np.array([list(atom_pos.value_in_unit_system(unit.md_unit_system)) for atom_pos in new_positions]), unit = unit.nanometers)

# Move old pdb
process = subprocess.run(['mv', args.name, f"{args.name[:-4]}_old.pdb"])

# Write topology and positions to pdb
with open(args.name, 'w') as f:
    app.PDBFile.writeFile(new_topology, new_positions_corrected, f, keepIds=True)
    

