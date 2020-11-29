## Manifest
* `01_openmm_integrator`- Contains scripts and PDB files used for setting up the first round of equilibration (using the OpenMM 7.4.2 Langevin integrator). Also contains equilibration script and output files.
* `02_openmmtools_integrator` - Contains scripts and PDB files used for setting up the second round of equilibration (using the OpenMMTools 0.20.0 BAOAB Langevin integrator), which was conducted because the BAOAB Langevin integrator has improved stability. Also contains equilibration script and output files, the latter of which was used to seed Folding@home simulations.
* `representative_glycan_structures` - Contains glycan structures from [Elisa Fadda] (https://www.maynoothuniversity.ie/people/elisa-fadda) / Aoife M. Harbison that were added to RBD:ACE2
