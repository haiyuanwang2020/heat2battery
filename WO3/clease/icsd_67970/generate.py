from ase.io import read
from clease.settings import ClusterExpansionSettings, Concentration
from ase.db import connect
from clease import NewStructures

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

# Generate and insert initial structures into database
ns = NewStructures(settings=settings, generation_number=0, struct_per_gen=50)
ns.generate_initial_pool()
