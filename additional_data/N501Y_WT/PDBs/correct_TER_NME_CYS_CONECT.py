import argparse

# Read filename
parser = argparse.ArgumentParser(description='remove CONECT records, correct NMEs and CYSs')
parser.add_argument('name', type=str, help='name of file for which to correct')
args = parser.parse_args()

# Read lines
with open(args.name, "r") as f:
    lines = f.readlines()

# Iterate through lines, copying them over to new list of lines
glycan_residue_names = ['UYB', '4YB', 'VMB', '2MA', '0YB', '0FA', '0LB']
new_lines = []
previous_res_id =  0
previous_res_name = ''
for line in lines:
    if 'CONECT' in line: # Skip CONECT lines
        continue
    if 'TER' not in line and 'END' not in line and 'REMARK' not in line:
        current_res_name = line[17:20]
        current_res_id = int(line[23:26])
        if current_res_id != previous_res_id:
            if previous_res_name in glycan_residue_names:
                new_lines.append("TER\n") # add TER if the previous residue was a glycan residue
            if current_res_name == "NME":
                new_lines = new_lines[:len(new_lines)-1] # remove the TER before the NME
            if previous_res_name == "NME": 
                new_lines.append("TER\n") # add TER after the NME and before starting the next residue
            previous_res_id = current_res_id
            previous_res_name = current_res_name
        if current_res_name == 'NME': # change C atom in NMEs to CH3
            atom = line[13:16]
            if atom == 'C  ':
                    line = line[:13] + 'CH3 ' + line[17:]
        if "RBD" in args.name:
            if current_res_name == 'CYS': # change CYS to CYX
                line = line[:17] + 'CYX' + line[20:]
    new_lines.append(line)

with open(args.name[:-4] + "_final.pdb", 'w') as f:
    f.writelines(new_lines)
