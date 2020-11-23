import argparse

# Read filename
parser = argparse.ArgumentParser(description='correct TERs')
parser.add_argument('name', type=str, help='name of file for which to correct')
args = parser.parse_args()

# Read lines
with open(args.name, "r") as f:
	lines = f.readlines()

# Iterate through lines, copying them over to new list of lines
new_lines = []
glycan_residue_names = ['UYB', '4YB', 'VMB', '2MA', '0YB', '0fA', '0LB']
CYS_to_keep = [261, 498] # List of CYS residues to keep (do not change these to CYX)
previous_res_id = 0
previous_res_name = ''
for line in lines:
    if 'TER' not in line and 'END' not in line:
        current_res_id = int(line[23:26])
        current_res_name = line[17:20]
        if current_res_id != previous_res_id:
            if previous_res_name in glycan_residue_names:
                new_lines.append("TER\n") # add TER if the previous residue was a glycan residue
            if current_res_name == "NME":
                new_lines = new_lines[:len(new_lines)-1] # remove the TER before the NME
            if previous_res_name == "NME": 
                new_lines.append("TER\n") # add TER after the NME and before starting the next residue
            previous_res_id = current_res_id
            previous_res_name = current_res_name
    new_lines.append(line)

with open(args.name[:-4] + "_corrected.pdb", 'w') as f:
    f.writelines(new_lines)
