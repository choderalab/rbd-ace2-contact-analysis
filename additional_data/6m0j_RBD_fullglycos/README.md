### Author: Alex Payne & William Glass

# System Preparation

## Structure Download
This pipeline uses Tristan Croll's refined structure of PDB:6m0j from the [CORONAVIRUS STRUCTURAL TASKFORCE DATABASE](https://github.com/thorn-lab/coronavirus_structural_task_force). `6m0j_refine_14.pdb` was downloaded from the [database](https://github.com/thorn-lab/coronavirus_structural_task_force/blob/master/pdb/surface_glycoprotein/SARS-CoV-2/6m0j/isolde/6m0j_refine_14.pdb). \

## Obtain starting structure

Location: `./PDBs`
`6m0j_refine_14.pdb` contains the RBD - ACE2 complex. On the RBD there are two GlcNAc residues at N343, as well as a core fucose, as modeled by Tristan Croll.

In PyMol, the RBD from `6m0j_refine_14`, as well as the nearby crystal waters, were extracted. The glycans were retained and used as a scaffold to align the FA2G2 glycan to. Additionally, the protein termini were capped with ACE and NME residues.

The core GlcNacs and Fucose in FA2G2 were aligned using the `pair_fit` tool in PyMol. Once aligned, the glycan was checked for any steric clashes with the protein side chains, of which there were none. The original GlcNac-(Fuc)-GlcNAc core was then removed from the structure leaving the RBD and aligned FA2G2 molecule.

This model was saved as `RBD_6m0j_refine_14_capped_N343glycosylated.pdb`.


## Cleaning the PDB file

All `TER` cards from `RBD_6m0j_refine_14_capped_N343glycosylated.pdb` were removed and `pdb_reatom` and `pdb_reres` (from [PDB Tools](http://www.bonvinlab.org/pdb-tools/)) used to renumber atoms and residues respectively. The `TER` cards need to be removed for this to work properly.

```bash
pdb_reatom RBD_6m0j_refine_14_capped_N343glycosylated.pdb > temp.pdb
pdb_reres -332 temp.pdb > temp2.pdb
mv temp2.pdb RBD_6m0j_refine_14_capped_N343glycosylated.pdb
```

The `TER` cards were then added at the end of the main protein chain, as well as between each glycan residue. Crystal water resides were moved to the end of the PDB file so that the final order of the PDB was: Protein -> Glycans -> Solvent.

The resulting file was saved as `RBD_6m0j_refine_14_capped_N343glycosylated_cleaned.pdb`.

## Running tleap

Location: `./run_tleap`

The PDB file must be pre-processed before using it in the AMBER `tleap` program. Namely, all glycosylated aspargine residues must be renamed from ASN -> NLN and cysteines involved in disulphide bridges need to be renamed from CYS -> CYX.

This was carried out on the `RBD_6m0j_refine_14_capped_N343glycosylated_cleaned.pdb` file and saved as `./run_tleap/RBD_6m0j_refine_14_capped_N343glycosylated_cleaned_LEAP_INPUT.pdb`.


### Bonding

Bonding between disulphide bridges needs to be specified in `tleap` as well as bonding between NLN residues and GlcNAc + the entire glycan complex. The details of bonds specified can be found in the `RBD_leap.in` file.

### Addition of solvent and ions

Within the `RBD_leap.in` file the `solvateBox` command solvates the system and the `addIonsRand` command replaces water molecules with Na+ and Cl- ions. 

The `numWaters` variable was determined by running the `RBD_tleap.in` script *without* adding ions first. The above method gave `numPositive` (i.e. Na+) as 33 and `numNegative` (i.e. Cl-) as 35.

The AMBER `tleap` program does not automatically work out the correct number of ions for a given system / system charge. In order for this to be determined we use the methodology from OpenMM:

```bash
from math import floor
numWaters = 12048
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
```

## Running tleap

The `tleap` command was run by: 

```bash
tleap -s -f RBD_tleap.in > RBD_tleap.out
```

This produced a fully solvated (with ions) system. Check `RBD_tleap.out` and `leap.log` for a full description of the output. The `.prmtop` and `.inpcrd` files can now be used for minimisation and equilibration.

## Equilibration

Once prepared the system was equilibrated here: `./equil`
