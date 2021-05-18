from gaussium.common import read_basis_set_file
from gaussium.common import read_mol_file
from gaussium.common import Symmetry
from gaussium.energy import Energy
from gaussium.factory import MoleculeFactory
from gaussium.geometryoptimization import NelderMead
import numpy as np
import time


def start(mol_file, basis_file, method, processors, symmetry=False, geometry_optimization=None):
    np.set_printoptions(linewidth=100000, threshold=np.inf)
    start_time = time.perf_counter()

    nuclei_list, electrons, multiplicity = read_mol_file(mol_file)
    energy_object = Energy(electrons, multiplicity, processors, method)

    print('\n*************************************************************************************************')
    print('\nA BASIC QUANTUM CHEMICAL PROGRAM IN PYTHON\n\n\n{}'.format([x.element for x in nuclei_list]))

    if geometry_optimization is not None:
        nelder_mead = NelderMead(basis_file, energy_object, nuclei_list)
        energy = nelder_mead.optimize()
    else:
        nuclei_list, point_group = MoleculeFactory(symmetry).create(nuclei_list)
        basis_set = read_basis_set_file(basis_file, nuclei_list)
        energy_object.symmetry_object = Symmetry(point_group, basis_set)
        energy = energy_object.calculate_energy(nuclei_list, basis_set)

    print('\n*************************************************************************************************')
    print('\nTIME TAKEN: ' + str(time.perf_counter() - start_time) + 's')
    print("\nWhat I cannot create I cannot understand - Richard Feynman\n")

    return energy
