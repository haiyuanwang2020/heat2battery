# cluster expansion of WO3 with hydrogens 

The ICSD 67970 is used here as a starting point with a stoichiometry of WO3H6. The primitive cell has 40 atoms (`POSCAR_primitive`).

1. Generate the initial configurations. Here we fix the W atoms, and allow the composition of O/H to vary between 0 (empty) or 1 (fully occupied).

`./generate.py`

2. Perform total energy calculations (here with MACE) with the generated configurations.

`./calc.py`

3. Fit the energies vs the correlation functions to obtain the ECIs.

`./eci.py`
