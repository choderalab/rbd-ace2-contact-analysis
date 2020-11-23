import argparse

# Read filename
parser = argparse.ArgumentParser(description='remove CONECT records, correct NMEs and CYSs')
parser.add_argument('name', type=str, help='name of file for which to correct')
args = parser.parse_args()

# Read lines
with open(args.name, "r") as f:
    lines = f.readlines()

# Iterate through lines, copying them over to new list of lines
CYS_to_keep = [261, 498] # List of CYS residues to keep (do not change these to CYX)
new_lines = []
for line in lines:
    if 'CONECT' in line: # Skip CONECT lines
        continue
    if 'TER' not in line and 'END' not in line and 'REMARK' not in line:
        current_res_name = line[17:20]
        current_res_id = int(line[23:26])
        if current_res_name == 'NME': # change C atom in NMEs to CH3
            atom = line[13:16]
            if atom == 'C  ':
                line = line[:13] + 'CH3 ' + line[17:]
        if current_res_name == 'CYS' and current_res_id not in CYS_to_keep: # change CYS to CYX
            line = line[:17] + 'CYX' + line[20:]
    new_lines.append(line)

with open(args.name[:-4] + "_final.pdb", 'w') as f:
    f.writelines(new_lines)
