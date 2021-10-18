### Author: William Glass

# Equilibration Protocol

The RBD:ACE2 N501Y system (`RBD_ACE2_complex_N501Y.inpcrd, RBD_ACE2_complex_N501Y.prmtop`) was equilibrated for 10 ns. Check `equil_NPT_10ns.py` for complete simulation details.

* Performed in the NPT ensemble using the OpenMMTools BAOAB Langevin Integrator and Monte Carlo Barostat.
* Use of hydrogen mass repartitioning (HMR). 
* Timestep of 4 fs.

Outputs are in `output/`

`.dcd` not included since it is large and unnecessary

Since tleap merges all chains into one and renumbers residues starting from 1, `output/equilibrated.pdb` was renumbered using `correct_equilibrated_rbd_ace2.py` and saved as `output/equilibrated.pdb`. The original, non-renumbered structure was moved to `output/equilibrated_old.pdb`.

Number of atoms: 200481
Time (s) per 10 ns: 20548.451
