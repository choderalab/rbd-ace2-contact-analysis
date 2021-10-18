### Author: Ivy Zhang

The RBD:ACE2 complex (fully glycosylated) was created using the prepared monomeric structures of the [RBD](https://github.com/choderalab/vir-antibody-structure-preparation/blob/main/systems/rbd_ace2_complex/WT/PDBs/3_refined_ace2_rbd_complex_RBD_ONLY_cleaned.pdb) and [ACE2](https://github.com/choderalab/vir-antibody-structure-preparation/blob/main/systems/rbd_ace2_complex/WT/PDBs/3_refined_ace2_rbd_complex_ACE2_ONLY_cleaned.pdb).

## Glycan structures

Glycosylation sites and patterns were determined by [Elisa Fadda](https://www.maynoothuniversity.ie/people/elisa-fadda). The recent [Science paper](https://science.sciencemag.org/content/early/2020/05/01/science.abb9983) by Watanabe *et al.* was used to propose the most likely glycosylation pattern:

Pre-equilibrated glycan structures are stored in the `system/representative_glycan_structures directory`.

# System Preparation

Since the monomeric systems had already been constructed, we first mutated N501Y in the RBD and then combined RBD N501Y and ACE2 to form the final RBD:ACE2 N501Y complex:

## Obtain starting structure

The WT RBD and ACE2 structure were copied from `../WT/PDBs` and saved here:

WT RBD location: `./PDBs/3_refined_ace2_rbd_complex_ACE2_ONLY_cleaned.pdb`

WT ACE2 location: `./PDBs/3_refined_ace2_rbd_complex_RBD_ONLY_cleaned.pdb`

Glycosylation patterns were kept the same as in the original structures, as a reminder these are:

**ACE2**
* N53 (~2% FA2, ~15% A3)
* N90 (~20% FA2[3/6]G1)
* N103 (~50% FA2)
* N322 (~1% FA2, ~60% other)
* N432 (~25% non-occupied, otherwise FA2G2 ok)
* N546 (15% A2, or A2[3/6]E1)
* N690 (~15% FA2)
* Only one of the three O-Glycan sites occupied: O730 but not in structure.

**RBD**
* N343 (FA2G2)

## Mutating N501Y in the RBD
The WT RBD structure was mutated to N501Y using PyMOL. The rotamer with the least number of clashes was chosen.

Location: `./PDBs/4_refined_ace2_rbd_complex_RBD_ONLY_N501Y.pdb`

## Cleaning the PDB file

The RBD N501Y structure was cleaned:`TER` cards were added between protein chains and glycan groups, NME residues were fixed (`TER` lines were moved to after the last NME atom), CYS residues were changed to CYX, and CONECT records were removed.

The final structure was saved in `./PDBs` as:

* `5_refined_ace2_rbd_complex_RBD_ONLY_N501Y_final.pdb`

Note, in the next section, the above RBD N501Y structure and the WT ACE2 structure (`./PDBs/3_refined_ace2_rbd_complex_ACE2_ONLY_cleaned.pdb`) are used.

## Running tleap

Location: `./run_tleap/`

Since all bonding between residues (such as disulphides and glycans) had already been specified in each monomeric system's `tleap` input file, these were used to construct the `rbd_ace2_complex_N501Y_tleap.in` input file. 

### Addition of solvent and ions

Within the `rbd_ace2_complex_N501Y_tleap.in` file the `solvateBox` command solvates the system and the `addIonsRand` command replaces water molecules with Na+ and Cl- ions. 

The AMBER `tleap` program does not automatically work out the correct number of ions for a given system / system charge. In order for this to be determined we use the methodology from OpenMM:

```python
from math import floor
numWaters = 61992
numPositive = 0
numNegative = 0 
totalCharge = -20
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

The `numWaters` variable was determined by running the `rbd_ace2_complex_N501Y_tleap.in` script *without* adding ions first. The above method gave `numPositive` (i.e. Na+) as 187 and `numNegative` (i.e. Cl-) as 167.


## Running tleap

The `tleap` command was run by: 

```bash
tleap -s -f rbd_ace2_complex_N501Y_tleap.in > rbd_ace2_complex_N501Y_tleap.out
```

This produced a fully solvated (with ions) system. Check `rbd_ace2_complex_N501Y_tleap.out` and `leap.log` for a full description of the output. The `.prmtop` and `.inpcrd` files could now be used for minimisation and equilibration. The generated system can be viewed here: `./run_tleap/RBD_ACE2_complex_N501Y.pdb`.

## Equilibration

Once prepared the system was equilibrated here: `./equil`
