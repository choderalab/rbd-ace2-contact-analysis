### Author: Alex Payne & William Glass

# System Preparation

## Obtain starting structure
Location: `./PDBs`
`RBD_6m0j_refine_14_capped.pdb` was taken from the adjacent `6m0j_RBD_fullglycos` directory, included in PDBs. A capped RBD structure from Tristan Croll's 6m0j_refine_14 structure. 

NAG and FUC residues were removed in PyMol and structure written out as `RBD_6m0j_refine_14_capped_noglycos.pdb`.

## Cleaning the PDB file

All `TER` cards from `RBD_6m0j_refine_14_capped_noglycos.pdb` were removed and `pdb_reatom` and `pdb_reres` (from [PDB Tools](http://www.bonvinlab.org/pdb-tools/)) used to renumber atoms and residues respectively. The `TER` cards need to be removed for this to work properly.

```bash
pdb_reatom RBD_6m0j_refine_14_capped_noglycos.pdb > temp.pdb
pdb_reres -332 temp.pdb > temp2.pdb
mv temp2.pdb RBD_6m0j_refine_14_capped_noglycos_renum.pdb
```

The `TER` cards were then added at the end of the main protein chain. 

The resulting file was saved as `RBD_6m0j_refine_14_capped_noglycos_cleaned.pdb`.

## Running tleap

Location: `./run_tleap`

The PDB file must be pre-processed before using it in the AMBER `tleap` program. Namely, all glycosylated aspargine residues must be renamed from ASN -> NLN and cysteines involved in disulphide bridges need to be renamed from CYS -> CYX.

This was carried out on the `RBD_6m0j_refine_14_capped_N343glycosylated_cleaned.pdb` file and saved as `RBD_6m0j_refine_14_capped_noglycos_cleaned_LEAP_INPUT.pdb` in `run_tleap`.

### Bonding

Bonding between disulphide bridges needs to be specified in `tleap`. The details of bonds specified can be found in the `RBD_leap.in` file.

### Addition of solvent and ions

Within the `RBD_leap.in` file the `solvateBox` command solvates the system and the `addIonsRand` command replaces water molecules with Na+ and Cl- ions. 

The `numWaters` variable was determined by running the `RBD_tleap.in` script *without* adding ions first. The method below gave `numPositive` (i.e. Na+) as 27 and `numNegative` (i.e. Cl-) as 29.

The AMBER `tleap` program does not automatically work out the correct number of ions for a given system / system charge. In order for this to be determined we use the methodology from OpenMM:

```python
from math import floor
numWaters = 10039
numPositive = 0
numNegative = 0 
totalCharge = 2
ionicStrength = 0.15

if totalCharge > 0:
    numNegative += totalCharge
else:
    numPositive -= totalCharge

numIons = (numWaters - numPositive - numNegative) * ionicStrength / (55.4)  # Pure water is about 55.4 molar (depending on temperature)
numPairs = int(floor(numIons + 0.5))
numPositive += numPairs
numNegative += numPairs
print(numPositive, numNegative)
```

## Running tleap

The `tleap` command was run by: 

```bash
tleap -s -f RBD_tleap.in > RBD_tleap.out
```

This produced a fully solvated (with ions) system. Check `RBD_tleap.out` and `leap.log` for a full description of the output. The `.prmtop` and `.inpcrd` files can now be used for minimisation and equilibration.

## Equilibration

Once prepared the system was equilibrated here: `./equil`