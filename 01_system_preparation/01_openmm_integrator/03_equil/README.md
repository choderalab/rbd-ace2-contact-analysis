### Author: William Glass

# Equilibration Protocol

The RBD:ACE2 system (`RBD_ACE2_complex_FullGlycos.inpcrd, RBD_ACE2_complex_FullGlycos.prmtop`) was equilibrated for 10 ns. Check `equil_NPT_10ns.py` for complete simulation details.

* Performed in the NPT ensemble using the Langevin Integrator and Monte Carlo Barostat.
* Use of hydrogen mass repartitioning (HMR). 
* Timestep of 4 fs.