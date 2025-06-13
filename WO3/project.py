from ase.io import read
from clease.settings import Concentration, CECrystal 

#atoms = read("EntryWithCollCode836.cif")

W1 = (0.25, 0.029, 0.031)
W2 = (0.25, 0.030, 0.532)
O1 = (0, 0, 0)   
O2 = (0, 0, 0.5)  
O3 = (0.25, 0.269, 0.027) 
O4 = (0.25, 0.278, 0.471) 
O5 = (0.25, 0.004, 0.262) 
O6 = (0.25, 0.015, 0.776) 
H1 = (0.6532,  0.0955,  0.09035)
H2 = (0.51063, 0.23503, 0.24974)
H3 = (0.08808, 0.1513,  0.1604)
H4 = (0.62565, 0.13425, 0.3552)
H5 = (0.6251,  0.3533,  0.39697)
H6 = (0.1223,  0.40262, 0.36525)
H7 = (0.15256, 0.42246, 0.09432)
H8 = (0.12943, 0.1245,  0.38714)
basis = [W1, W2, O1, O2, O3, O4, O5, O6, H1, H2, H3, H4, H5, H6, H7, H8]

# specify the concentration ranges of species: e.g. WO3Hx
conc = Concentration(basis_elements=[['W'], ['W'], ['O'], ['O'], ['O'], ['H'], ['H'], ['H'], ['H'], ['H'], ['H'], ['H'], ['H']],
                     grouped_basis=[[0, 1],[2, 3, 4, 5, 6, 7], [8, 9, 10, 11, 12, 13, 14, 15]]
                    )

conc.set_conc_ranges(ranges=[[(1, 1)], [(3, 3)], [(0, 1)]])

'''
conc.set_conc_formula_unit(
   # formulas=["W<x>O<y>H<z>"],
    formulas=["W<x>", "O<y>", "H<z>"],
    variable_range={"x": (1, 1), "y": (3, 3), "z": (0., 1.)}
   # variable_range={
   #     "x": (2, 2),
   #     "y": (6, 6),
   #     "z": (0, 8)
   # }
)
'''

settings = CECrystal(
    basis=basis,
    cell=[7.34, 7.57, 7.75, 90, 90, 90],
    spacegroup=62,
    concentration=conc,
    size=[2, 2, 2],            # controls how big your supercell is
    max_cluster_dia=[6.0, 4.5, 4.5],
    db_name="wo3hx.db"         # your database for storing structures/energies
)

# generating initial structures
generator = NewStructures(settings)
generator.generate_initial_pool(n_structures=20,
                                max_per_conc=5   # At most 5 structures for each concentration of H
                                )

setting.save('WO3Hx.json')


