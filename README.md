# RBD ACE2 contact analysis
Code and workflow for running MD simulations on Folding@home to analyze RBD:ACE2 contacts

Publication: 

## License
* This software is licensed under the [MIT license](https://opensource.org/licenses/MIT) - a copy of this license is provided as `SOFTWARE_LICENSE`
* The data in this repository is made available under the Creative Commons [CC0 (“No Rights Reserved”) License](https://creativecommons.org/share-your-work/public-domain/cc0/) - a copy of this license is provided as `DATA_LICENSE`

## Manifest

* `01_system_preparation` - Contains all scripts and PDBs before and after structure preparation, system parameterization, and equilibration. Also contains the representative glycan structures added to RBD:ACE2
* `02_analysis` - Contains the script used to analyze RBD:ACE2 contacts in the trajectories

## Contributors

* Ivy Zhang
* William G. Glass
* Tristan Croll
* Aoife M. Harbison
* Elisa Fadda
* John D. Chodera

# Simulation details

**Structure preparation**
The RBD:hACE2 complex was constructed from individual RBD (PDB: [6m0j](https://www.rcsb.org/structure/6M0J), Chain E) and hACE2 (PDB: [1r42](https://www.rcsb.org/structure/1R42), Chain A) monomers aligned to the full RBD:hACE2 structure (PDB: 6m0j). The 1r42 structure was used for hACE2 because 1) 1r42 is higher resolution (2.20 Å, whereas 6m0j is 2.45 Å) and 2) the electron density map of 1r42 clearly reveals the NAG orientation at each glycosylated asparagine residue, providing a reliable building block on which to construct more complex glycan structures. While the NAGs for each glycan were present in 1r42, we constructed more complex glycans at each NAG due to the important role of glycans in mediating RBD:hACE2 binding, as shown by [Casalino et al.](https://doi.org/10.1021/acscentsci.0c01056). 

In order to start from the most reliable structural models, we obtained 6m0j and 1r42 from the [Coronavirus Structural Taskforce (CST) database](https://github.com/thorn-lab/coronavirus_structural_task_force), which contains refined structural models based on careful examination of the electron density. In the [RBD of the refined 6m0j structure](https://github.com/thorn-lab/coronavirus_structural_task_force/blob/master/pdb/surface_glycoprotein/SARS-CoV-2/6m0j) amino acid rotamers and peptide bonds were flipped to increase Ramanchandran favourability, decrease rotamer outliers, and reduce clashes. A more detailed summary of the 6m0j refinement details is available at: https://github.com/thorn-lab/coronavirus_structural_task_force/blob/master/pdb/surface_glycoprotein/SARS-CoV-2/6m0j/isolde/notes.txt. The 1r42 refined structure differs from the PDB-deposited structure in that it includes the missing C-terminal domain of hACE2 (copied from the 6m17 PDB structure). A more detailed summary of the 1r42 refinement details is available [here](https://github.com/thorn-lab/coronavirus_structural_task_force/blob/master/pdb/surface_glycoprotein/SARS-CoV-2/6m0j/isolde/notes.txt).

The resulting RBD and hACE2 monomers were then aligned in [PyMOL 2.3.2](http://pymol.org) to the [CST 6m0j structure](https://github.com/thorn-lab/coronavirus_structural_task_force/tree/master/pdb/surface_glycoprotein/SARS-CoV-2/6m0j) to create an initial RBD:hACE2 complex. The overall RMSD was 0.426 Å and the interface RMSD was 0.405 Å, where RMSD was computed for all atoms and the interface residues were defined as all residues within 4 Å of the other binding partner.

Next, the full glycosylation patterns for hACE2 and RBD glycans were determined from [Shajahan et al.](http://dx.doi.org/10.1093/glycob/cwaa101) and [Watanabe et al](https://doi.org/10.1126/science.abb9983). For the constructed RBD:hACE2 complex, these included sites: N53, N90, N103, N322, N432, N546, and N690 on hACE2 and N343 on the RBD. The glycan structures used for each site (FA2, FA26G1, FA2, FA2, FA2G2, A2, FA2, FA2G2, respectively) correspond to the most stable conformers obtained from multi microsecond MD simulations of cumulative sampling ([Harbison et al.](https://pubmed.ncbi.nlm.nih.gov/30325416/)). Base NAG residues of each glycan structure were aligned to the corresponding NAG stub in the RBD:hACE2 model in PyMOL 2.3.2 and any resulting clashes were refined in [ISOLDE](https://isolde.cimr.cam.ac.uk/).

**System solvation and parametrization**
The refined glycosylated RBD:hACE2 complex was prepared for simulation using the [AmberTools17](https://ambermd.org/AmberTools.php) tleap suite. All relevant disulfide bridges were specified as well as covalent connectivity within each glycan structure. The glycosylated protein was parameterized with the Amber ff14SB ([Maier et al., 2015](https://doi.org/10.1021/acs.jctc.5b00255)) and GLYCAM_06j-1 ([Kirschner et al., 2008](https://doi.org/10.1002/jcc.20820)) force fields. The system was solvated using the TIP3P rigid water model ([Jorgensen et al., 1983](https://doi.org/10.1063/1.445869)) in a cubic box with 1.5 nm solvent padding on all sides. The solvated system was then minimally neutralized with 0.15 M NaCl using the Li/Merz ion parameters of monovalent ions for the TIP3P water model (12-6 normal usage set) (Li et al., 2015 PMID: 26574374). Full details and tleap scripts can be found at: https://github.com/choderalab/rbm-ace2-contact-analysis.

**System equilibration**
The system was energy-minimized with an energy tolerance of 10 kJ mol−1and equilibrated using the OpenMM 7.4.2 Langevin integrator for 300 ns in the NPT (p=1 atm, T = 310 K) ensemble with a timestep of 4.0 femtoseconds, a collision rate of 1.0 picoseconds-1, and a constraint tolerance of 1 ✕ 10−5. Hydrogen atom masses were set to 4.0 amu by transferring mass from connected heavy atoms, bonds to hydrogen were constrained, and center of mass motion was not removed. Pressure was controlled by a molecular-scaling Monte Carlo barostat with an update interval of 25 steps. Non-bonded interactions were treated with the Particle Mesh Ewald method ([Darden et al., 1993](https://doi.org/10.1063/1.464397)) using a real-space cutoff of 1.0 nm and the OpenMM default relative error tolerance of 0.0005, with grid spacing selected automatically. For improved stability, the structure was then equilibrated using the [OpenMMTools 0.20.0](http://openmmtools.org) BAOAB Langevin integrator ([Leimkuhler and Matthews, 2013](https://aip.scitation.org/doi/10.1063/1.4802990)) for 10 ns using all of the same simulation parameters described above. This simulation was subsequently packaged to seed for production simulation on Folding@home ([Shirts and Pande, 2000](https://science.sciencemag.org/content/290/5498/1903.full), [Zimmerman et al., 2020](https://doi.org/10.1101/2020.06.27.175430)). Default parameters were used unless noted otherwise. Further details of the equilibration protocol are available at: https://github.com/choderalab/rbm-ace2-contact-analysis

**Folding@home simulation**
The equilibrated structure was then used to initiate parallel distributed MD simulations on Folding@home ([Shirts and Pande, 2000](https://science.sciencemag.org/content/290/5498/1903.full), [Zimmerman et al., 2020](https://doi.org/10.1101/2020.06.27.175430)). Simulations were run with OpenMM 7.4.2 (Folding@home core22 0.0.13). Production simulations used the same Langevin integrator as the NpT equilibration described above. In total, 2000 independent MD simulations were generated on Folding@home. Conformational snapshots (frames) were stored at an interval of 0.5 ns/frame for subsequent analysis. The resulting final dataset contained 2000 trajectories, 173.8 µs of aggregate simulation time, and 34761 frames. This amount of simulation time corresponds to ~12.9 GPU-years on an NVIDIA GeForce GTX 1080Ti. This trajectory dataset with solvent is available at the MolSSI COVID-19 Molecular Structure and Therapeutics Hub: https://covid.molssi.org//org-contributions/#folding--home (**TO BE UPDATED**)
