from ase.db import connect
from mace.calculators import MACECalculator
from clease.tools import update_db
from ase.calculators.singlepoint import SinglePointCalculator
import numpy as np

models = '/gpfs/home/acad/ucl-modl/wchen/mace_models/mace-omat-0-medium.model'
db_name = 'wo3.db'

calc = MACECalculator(model_paths=models, device='cuda')
db = connect(db_name)

for row in db.select(converged=False):
    atoms = row.toatoms()
    atoms_no_X = atoms.copy()
    non_vacancy_mask = np.array([s != 'X' for s in atoms.get_chemical_symbols()])
    atoms_no_X = atoms[non_vacancy_mask]

    atoms_no_X.calc = calc
    energy = atoms_no_X.get_potential_energy()

    sp_calc = SinglePointCalculator(atoms, energy=energy)
    atoms.calc = sp_calc
    update_db(uid_initial=row.id, final_struct=atoms, db_name=db_name)
