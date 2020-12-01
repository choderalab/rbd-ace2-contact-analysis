### Author: William Glass

# Equilibration Protocol

The RBD:ACE2 system (`../02_run_tleap/RBD_ACE2_complex_FullGlycos.inpcrd`, `../02_run_tleap/RBD_ACE2_complex_FullGlycos.prmtop`) was equilibrated for 10 ns. Check `equil_NPT_10ns.py` for complete simulation details. A brief outline is given below:

* Performed in the NPT ensemble using the `OpenMM` Langevin Integrator and Monte Carlo Barostat.
* Use of hydrogen mass repartitioning (HMR). 
* Timestep of 4 fs.

The equilibrated system (PDB and XML files) is located in `./output`.
