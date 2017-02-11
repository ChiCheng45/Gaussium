import sys, os
sys.path.insert(1, os.path.dirname(os.path.realpath(__file__)) + '/../../')

from src.main.common import read_basis_set_file
from src.main.common import read_mol_file
from src.main.common import Symmetry
from src.main.common import Energy
from src.main.objects import Molecule
from src.main.objects import PointGroup
from src.main.factory import MoleculeFactory
import numpy as np
import time


def menu():
    # start('HeH+.mol', 'STO-3G.gbs', 'RHF', 4)  # -2.84183608212 a.u.
    # start('HeH+.mol', '6-311+GPP.gbs', 'RHF', 4)  # -2.92922773384 a.u.
    # start('C2H4.mol', '3-21G.gbs', 'RHF', 4)  # -77.600460844 a.u. 30.747198048700866s
    # start('O2.mol', 'STO-3G.gbs', 'UHF', 4, symmetry=True)  # -147.634028141 a.u.
    # start('O2.mol', 'STO-3G.gbs', 'GUHF', 4)  # -147.634028141 a.u.
    # start('CO.mol', 'STO-3G.gbs', 'MP2', 4)  # -111.354512528 a.u.
    # start('H2O.mol', 'STO-3G.gbs', 'RHF', 4, symmetry=True)
    start('C2H4.mol', '3-21G.gbs', 'RHF', 4, symmetry=True)  # -77.600460844 a.u. 19.0269839632222s
    # start('H2O.mol', 'STO-3G.gbs', 'CIS', 4)  # 0.2872554996 a.u. 0.3564617587 a.u.

    # only worth doing DFT calculations on atoms at the moment
    # start('He.mol', 'STO-3G.gbs', ('DFT', 'S', ''), 4)  # -2.657311972 a.u.
    # start('He.mol', 'STO-3G.gbs', ('DFT', 'S', 'VWN3'), 4)  # -2.80983127318 a.u.
    # start('Li-.mol', 'STO-3G.gbs', ('DFT', 'S', 'VWN3'), 4)  # -7.223380048745456 a.u.

    # start('H2O.mol', 'STO-3G.gbs', 'CCSD', 4)  # -0.0706800939192 a.u.
    # start('CH4.mol', 'STO-3G.gbs', 'CCSD', 4)  # -0.078469894846414 a.u.
    # start('H2O.mol', 'STO-3G.gbs', 'CCSD(T)', 4)  # -9.98772699528e-05 a.u.


def start(mol_file, basis_file, method, processors, symmetry=False, geometry_optimization=None):
    np.set_printoptions(linewidth=100000, threshold=np.inf)
    start_time = time.clock()
    energy = 0.0

    nuclei_array, electrons, multiplicity = read_mol_file(mol_file)

    if symmetry:
        molecule = MoleculeFactory().point_group(nuclei_array)
    else:
        molecule = Molecule(nuclei_array, PointGroup([], [], [], [], 'C_{1}'))

    basis_set = read_basis_set_file(basis_file, molecule.nuclei_array)
    symmetry_object = Symmetry(molecule.point_group, basis_set)

    print('\n*************************************************************************************************')
    print('\nA BASIC QUANTUM CHEMICAL PROGRAM IN PYTHON\n\n\n{}'.format([x.element for x in nuclei_array]))

    if geometry_optimization is not None:
        energy_object = Energy(electrons, multiplicity, processors, method)
    else:
        energy = Energy(
            electrons, multiplicity, processors, method
        ).calculate_energy(molecule.nuclei_array, basis_set, symmetry_object)

    print('\n*************************************************************************************************')
    print('\nTIME TAKEN: ' + str(time.clock() - start_time) + 's')
    print("\nWhat I cannot create I cannot understand - Richard Feynman\n")

    return energy

if __name__ == "__main__":
    menu()
