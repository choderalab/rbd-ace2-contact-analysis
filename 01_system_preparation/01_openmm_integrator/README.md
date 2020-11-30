### Author: William Glass

The RBD:ACE2 complex (fully glycosylated) was created using the prepared monomeric structures of the RBD and ACE2 (see `../01_PDBs/3_refined_ace2_rbd_RBD_ONLY_cleaned.pdb` and `../01_PDBs/3_refined_ace2_rbd_ACE2_ONLY_cleaned.pdb`).

## Glycan structures

Glycosylation sites and patterns were determined by [Elisa Fadda](https://www.maynoothuniversity.ie/people/elisa-fadda). The [recent work](https://science.sciencemag.org/content/early/2020/05/01/science.abb9983) by Watanabe *et al.* and Shajahan *et al.* was used to propose the most likely glycosylation pattern in the RBD and ACE2 respectively:

Pre-equilibrated glycan structures, prepared by Aoife M. Harbison, are stored in the `../../structures-in-progress/representative_glycan_structures` directory.

Glycosylation patterns were determined as:

**ACE2**
* N53 (~2% FA2, ~15% A3)
* N90 (~20% FA2[3/6]G1)
* N103 (~50% FA2)
* N322 (~1% FA2, ~60% other)
* N432 (~25% non-occupied, otherwise FA2G2 ok)
* N546 (15% A2, or A2[3/6]E1)
* N690 (~15% FA2)

**RBD**
* N343 (FA2G2)
# System Preparation

## Addition of glycans to ACE2 and RBD
### ACE2
The refined ACE2 model was obtained from a refined `1r42` crystal structure prepared by Tristan Croll. In PyMol, protein termini were capped with ACE and NME residues. All GlcNAc / NAG residues were present and therefore kept as scaffold from which to place the complex glycan structures. Using the `Pair Fitting` wizard in PyMol, the base GlcNac / NAG residues in the appropriate glycan structures (`../../structures-in-progress/representative_glycan_structures`) were superimposed with the base GlcNaC / NAG residues in the refined `1r42` structure at each glycosylation site. Once aligned, the base GlcNac / NAG residue in the refined `1r42` structure was removed to leave only the complex glycan at each site.

Once all complex glycans were added the structure was checked for steric clashes. All sites had no major clashes, with the exception of the glycan at N53. Here, the base Fucose clashed with side chains of the protein and as a result it was manually rotated to remove severe clashes.

The `1r42` structure was used for hACE2 because:
* `1r42` is higher resolution (2.20 Å, whereas `6m0j` is 2.45 Å)
* The electron density map of `1r42` clearly reveals the NAG orientation at each glycosylated asparagine residue, providing a reliable building block on which to construct more complex glycan structures.

### RBD
The same procedure as detailed above was followed for the RBD at N343 with the FA2G2 glycan. No major clashes were observed between the glycan and protein.
## Obtaining ACE:RBD complex structure

Location: `./01_PDBs/0_RBD_ACE2_complex_from_1r42ref_6m0j_FyllGlycos_alignto6m0jRefined_NOCAPS.pdb`

`./01_PDBs/0_RBD_ACE2_complex_from_1r42ref_6m0j_FyllGlycos_alignto6m0jRefined_NOCAPS.pdb` was constructed by aligning the glycosylated `6m0j` structure RBD (above) and the refined `1r42` structure ACE2 (above) with the `6m0j` [CST structure](https://github.com/thorn-lab/coronavirus_structural_task_force/tree/master/pdb/surface_glycoprotein/SARS-CoV-2/6m0j) RBD:ACE2 complex. The resulting aligned proteins were then written out to file in PyMol.
## ISOLDE

After visual inspection of the ACE2:RBD interface it was clear there were severe steric clashes between ACE2 and RBD protein chains. In order to refine the structure ISOLDE was used. For an in-depth discussion of steps taken within the ISOLDE package check `./01_PDBs/README.md`.

## Cleaning the PDB file

The refined complex structure (`./PDBs/2_refined_ace2_rbd_complex.pdb`) was split into the ACE2 and RBD parts. This was done so that preparation with `tleap` (see later) could be performed easily. In each file the order of residues was changed so that each had the following order: protein -> glycans -> ions (if present). Then, `TER` cards were added between protein chains and glycan groups. Atoms were renumbered using `pdb_reatom` (from [PDB Tools](http://www.bonvinlab.org/pdb-tools/)).

The final structures are saved in `./01_PDBs` as:

* `3_refined_ace2_rbd_complex_ACE2_ONLY_cleaned.pdb`
* `3_refined_ace2_rbd_complex_RBD_ONLY_cleaned.pdb` 


Both `3_refined_ace2_rbd_complex_ACE2_ONLY_cleaned.pdb` and `3_refined_ace2_rbd_complex_RBD_ONLY_cleaned.pdb` were then copied to `../02_run_tleap`.

## Running tleap

Location: `./02_run_tleap/`

Bonding between residues (such as disulphides and glycans etc) needed to be specified for both ACE2 and the RBD. The bonding patterns are explicitly specified in the `rbd_ace2_complex_tleap.in` input file. 
### Addition of solvent and ions

Within the `rbd_ace2_complex_tleap.in` file the `solvateBox` command solvates the system and the `addIonsRand` command replaces water molecules with Na+ and Cl- ions. 

The AMBER `tleap` program does not automatically work out the correct number of ions for a given system / system charge. In order for this to be determined the methodology from OpenMM was used, an example of which is outlined below:

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

The `numWaters` variable was determined by running the `rbd_ace2_complex_tleap.in` script *without* adding ions first. The above method gave values for `numPositive` (i.e. Na+) and `numNegative` (i.e. Cl-) that were then used in `rbd_ace2_complex_tleap.in` to add the required number of ions.
## Running tleap

The `tleap` command was run by: 

```bash
tleap -s -f rbd_ace2_complex_tleap.in
```

This produced a fully solvated (with ions) system. The `.prmtop` and `.inpcrd` files could now be used for minimisation and equilibration.
## Equilibration

Once prepared the system was equilibrated here: `./03_equil`