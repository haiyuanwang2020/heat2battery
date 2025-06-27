from clease.settings import ClusterExpansionSettings, Concentration
from ase.db import connect
from ase.io import read
import numpy as np
import json
from clease.montecarlo.sgc_montecarlo import SGCMonteCarlo
from clease.calculator import attach_calculator
from clease.montecarlo.observers import Snapshot

primitive_atoms = read('./POSCAR_primitive')

# Define concentration ranges explicitly
conc = Concentration(basis_elements=[
    ('W',),           # Fixed tungsten sublattice (only W)
    ('O', 'X'),       # Oxygen sublattice (occupied or vacancy)
    ('H', 'X')        # Hydrogen sublattice (occupied or vacancy)
])

conc.set_conc_ranges(ranges=[[(1, 1)], [(0, 1), (0, 1)], [(0, 1), (0, 1)]])

conc.get_random_concentration()

settings = ClusterExpansionSettings(
    prim=primitive_atoms,
    concentration=conc,
    supercell_factor=3,  
    db_name="wo3.db",
    max_cluster_dia=[6.0, 3.5],
    basis_func_type="polynomial",
)

# Load the ECIs
with open('eci.json', 'r') as f:
    eci = json.load(f)

# Create a large supercell for MC simulations
atoms = settings.atoms.copy()*(4,4,4)
atoms = attach_calculator(settings, atoms=atoms, eci=eci)

snap = Snapshot(atoms, fname='snap.traj')
mc = SGCMonteCarlo(atoms, symbols=['O','H','X'], temp=300)
mc.attach(snap, interval=100)
# Specify the chemical potential for the two species 
# (applied to the respective point correlation function)
mc.run(steps=1000, chem_pot={'c1_0':-0.5, 'c1_1':-0.3})

# Get results
thermo = mc.get_thermodynamic_quantities()
print(thermo)
