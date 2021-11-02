### Author: William Glass and Alex Payne

# Equilibration Protocol

The RBD system (`RBD_6m0j_refine_14_NoGlycos.inpcrd, RBD_6m0j_refine_14_NoGlycos.prmtop`) was equilibrated for 10 ns. Check `equil_NPT_10ns.py` for complete simulation details.

* Performed in the NPT ensemble using the OpenMMTools BAOAB Langevin Integrator and Monte Carlo Barostat.
* Use of hydrogen mass repartitioning (HMR). 
* Timestep of 4 fs.

Outputs (in `output/`):

`.dcd` final was not added to this directory as it is large and unnecessary

Since `tleap` merges all chains into one and renumbers residues starting from 1, `outp
ut/equilibrated.pdb` was renumbered using `correct_equilibrated_rbd_noglycos.py` and
 saved as `output/equilibrated.pdb`. The original, non-renumbered structure was moved
to `output/equilibrated_old.pdb`.

* Number of atoms: 39868
* Time (s) per 10 ns: 5285.201
