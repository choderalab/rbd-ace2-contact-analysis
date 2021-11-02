import MDAnalysis as mda
import argparse

parser = argparse.ArgumentParser(description='run mda to fix glycan vis')
parser.add_argument('pdb', type=str, help='path to input file')
args = parser.parse_args()

ref = mda.Universe(args.pdb)
protein = ref.select_atoms("all")
protein.write(args.pdb[:-4] + "_mda.pdb")
