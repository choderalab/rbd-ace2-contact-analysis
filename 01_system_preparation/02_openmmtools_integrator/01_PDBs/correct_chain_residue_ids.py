from simtk.openmm import app
import argparse

# Read filename
parser = argparse.ArgumentParser(description='correct TERs')
parser.add_argument('name', type=str, help='name of file for which to correct TERs')
parser.add_argument('--split', default=False, action='store_true')
args = parser.parse_args()
pdb = app.PDBFile(args.name)

# Rename chains
for chain in pdb.topology.chains():
    if chain.index == 0:
        chain.id = 'E'
        chain.index = 2
    elif chain.index >= 1 and chain.index <=10:
        chain.id = 'X'
        chain.index = 4
    elif chain.index == 11:
        chain.id = 'A'
        chain.index = 0
    elif chain.index >= 12 and chain.index <= 69:
        chain.id = 'X'
        chain.index = 3
    elif chain.index == 70:
        chain.id = 'A'
        chain.index = 1

# Get first residues of each chain
d_first_residues = {}
for residue in pdb.topology.residues():
    if "A" not in d_first_residues and residue.chain.id == "A":
        d_first_residues["A"] = residue.id
    elif "E" not in d_first_residues and residue.chain.id == 'E':
        d_first_residues["E"] = residue.id
    elif "X" not in d_first_residues and residue.chain.id == 'X':
        d_first_residues["X"] = residue.id
        
# Rename residues
current_chain = ""
d_current_start = {"A": 18, "E": 332, "X": 729}
for residue in pdb.topology.residues():
    if residue.chain.id != current_chain: # If the residue's chain does not match the current_chain 
        current_chain = residue.chain.id 
    if residue.id != d_first_residues[current_chain]: # If it's the first residue of the chain
        d_current_start[current_chain] += 1
    # Update the residue id
    old_id = residue.id
    residue.id = str(d_current_start[current_chain])

if not args.split:
	app.PDBFile.writeFile(pdb.topology, pdb.positions, open(f"{args.name[:-4]}_renumbered.pdb", "w"), keepIds=True)
else:
	# Get lists of chains to delete for each 
	rbd_to_delete = []
	ace2_to_delete = []
	for chain in pdb.topology.chains():
	    if chain.index in [0, 1, 3]:
	        rbd_to_delete.append(chain)
	    else:
	        ace2_to_delete.append(chain)

	modeller_rbd = app.Modeller(pdb.topology, pdb.positions)
	modeller_rbd.delete(rbd_to_delete)

	app.PDBFile.writeFile(modeller_rbd.topology, modeller_rbd.positions, open(f"{args.name[:-4]}_renumbered_RBD.pdb", "w"), keepIds=True)

	modeller_ace2 = app.Modeller(pdb.topology, pdb.positions)
	modeller_ace2.delete(ace2_to_delete)

	app.PDBFile.writeFile(modeller_ace2.topology, modeller_ace2.positions, open(f"{args.name[:-4]}_renumbered_ACE2.pdb", "w"), keepIds=True)