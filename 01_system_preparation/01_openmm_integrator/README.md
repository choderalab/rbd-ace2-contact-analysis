### Author: William Glass

The RBD:ACE2 complex (fully glycosylated) was created using the prepared monomeric structures of the RBD and ACE2 (see `../rbd_systems/6m0j_RBD_fullglycos` and `../ACE2_only`).

## Glycan structures

Glycosylation sites and patterns were determined by [Elisa Fadda](https://www.maynoothuniversity.ie/people/elisa-fadda). The recent [Science paper](https://science.sciencemag.org/content/early/2020/05/01/science.abb9983) by Watanabe *et al.* was used to propose the most likely glycosylation pattern:

Pre-equilibrated glycan structures are stored in the `structures-in-progress/representative_glycan_structures` directory. 

# System Preparation

Since the monomeric systems had already been constructed these were combined to form the final RBD:ACE2 complex:

## Obtain starting structure

Location: `./PDBs/0_RBD_ACE2_complex_from_1r42ref_6m0j_FyllGlycos_alignto6m0jRefined_NOCAPS.pdb`

`0_RBD_ACE2_complex_from_1r42ref_6m0j_FyllGlycos_alignto6m0jRefined_NOCAPS.pdb` was constructed by aligning the RBD system (`../rbd_systems/run_tleap/6m0j_RBD_fullglycosRBD_6m0j_refine_14_capped_N343glycosylated_cleaned_LEAP_INPUT.pdb`) and the ACE2 system (`../ACE2_only/run_tleap/1r42_refine_2_capped_FullGlycos_cleaned_LEAP_INPUT.pdb`) with the ISOLDE refined `6m0j` PDB RBD:ACE2 complex structure. The resulting aligned proteins were then written out to file in PyMol.

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

## ISOLDE

After visual inspection of the ACE2:RBD interface it was clear there were severe steric clashes between protein chains. In order to refine the structure the ISOLDE tool was used. For an in-depth discussion of steps taken within the ISOLDE package check `./PDBs/README.md`.

## Cleaning the PDB file

The refined complex structure (`./PDBs/2_refined_ace2_rbd_complex.pdb`) was split into the ACE2 and RBD parts. This was done so that preparation with `tleap` (see later) could be performed easily. In each file the order of residues was changed so that each had the following order: protein -> glycans -> ions (if present). Then, `TER` cards were added between protein chains and glycan groups. Atoms were renumbered using `pdb_reatom` (from [PDB Tools](http://www.bonvinlab.org/pdb-tools/)).

The final structures are saved in `./PDBs` as:

* `3_refined_ace2_rbd_complex_ACE2_ONLY_cleaned.pdb`
* `3_refined_ace2_rbd_complex_RBD_ONLY_cleaned.pdb` 


Both `3_refined_ace2_rbd_complex_ACE2_ONLY_cleaned.pdb` and `3_refined_ace2_rbd_complex_RBD_ONLY_cleaned.pdb` were then copied to `../run_tleap`.

## Running tleap

Location: `./run_tleap/`

Since all bonding between residues (such as disulphides and glycans) had already been specified in each monomeric system's `tleap` input file, these were used to construct the `rbd_ace2_complex_tleap.in` input file. 

### Addition of solvent and ions

Within the `rbd_ace2_complex_tleap.in` file the `solvateBox` command solvates the system and the `addIonsRand` command replaces water molecules with Na+ and Cl- ions. 

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

The `numWaters` variable was determined by running the `rbd_ace2_complex_tleap.in` script *without* adding ions first. The above method gave `numPositive` (i.e. Na+) as 187 and `numNegative` (i.e. Cl-) as 167.


## Running tleap

The `tleap` command was run by: 

```bash
tleap -s -f rbd_ace2_complex_tleap.in > rbd_ace2_complex_tleap.out
```

This produced a fully solvated (with ions) system. Check `rbd_ace2_complex_tleap.out` and `leap.log` for a full description of the output. The `.prmtop` and `.inpcrd` files could now be used for minimisation and equilibration. The generated system can be viewed here: `./run_tleap/RBD_ACE2_complex_FullGlycos.pdb`.

## Equilibration

Once prepared the system was equilibrated here: `./equil`