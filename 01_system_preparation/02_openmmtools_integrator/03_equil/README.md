### Author: Ivy Zhang

## Equilibration Protocol

The 4 RBD:ACE2 WT and mutant systems were equilibrated for 10 ns. Each of the input/output files have been renamed to `0-3` so that the simulations could be run using a batch job, see table below for the mappings. Check `equil_NPT_10ns.py` for complete simulation details.

- Performed in the NPT ensemble using the Openmmtools Langevin Integrator and Monte Carlo Barostat.
- Use of hydrogen mass repartitioning (HMR).
- Timestep of 4 fs.
Outputs are in `output/`

Note: `.dcd` not included since it is large and unnecessary

## Info for Fah

| Simulation number |  System name                        | Number of atoms | Time (s) per 10 ns | 
| ----------------- | ----------------------------------- | --------------- | ------------------ |
| 0                 | RBD:ACE2 N439K CLONE2 GEN58         | 207348          | 18727.216          |
| 1                 | RBD:ACE2 K417V CLONE2 GEN58         | 207342          | 23365.133          |
| 2                 | RBD:ACE2 N439K / K417V CLONE2 GEN58 | 207340          | 23091.160          |
| 3                 | RBD:ACE2 WT CLONE2 GEN58            | 207350          | 23425.071          |
