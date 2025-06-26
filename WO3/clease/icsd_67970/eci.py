from clease.settings import ClusterExpansionSettings, Concentration
from clease import Evaluate
from clease.corr_func import CorrFunction
import clease.plot_post_process as pp
from ase.db import connect
from ase.io import read
import numpy as np
import matplotlib.pyplot as plt

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

# Calculate correlation functions of the clusters
cf = CorrFunction(settings)
# Reconfigure the correlation function values of the entries in DB.
cf.reconfigure_db_entries(verbose=True)

# Fit the ECIs 
select_cond=[('converged','=','True')]
evl=Evaluate(settings,select_cond=select_cond,scoring_scheme='loocv', max_cluster_size=4)
# Fit with the optimal value of the hyperparameter alpha
alphas = evl.alpha_CV()
alpha = alphas[0][np.argmin(alphas[1])]
evl.set_fitting_scheme(fitting_scheme='l2', alpha=alpha)
evl.fit()
# Save ECIs to a json file
evl.save_eci(fname="eci.json")

# Plot E_pred vs E_dft(mlip)
fig = pp.plot_fit(evl, interactive=True)
plt.show()

