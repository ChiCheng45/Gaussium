from src.main.common import read_basis_set_file
from src.main.common import read_mol_file
from src.main.common import coulomb_matrix
from src.main.common import Symmetry
from src.main.objects import Molecule
from src.main.objects import MoleculeFactory
from src.main.hartreefock import RestrictedHF
from src.main.hartreefock import DODSUnrestricted
from src.main.hartreefock import ConstrainedUnrestricted
from src.main.hartreefock import BlockedHartreeFock
from src.main.moellerplesset import MoellerPlesset
import numpy as np
import time


def menu():
    # start('HeH+.mol', 'STO-3G.gbs', 'RHF')  # -2.84183608212 a.u.
    # start('HeH+.mol', '6-311+GPP.gbs', 'RHF')  # -2.92922773384 a.u.
    # start('C2H4.mol', '3-21G.gbs', 'RHF')  # -77.600460844 a.u. 30.747198048700866s
    # start('O2.mol', 'STO-3G.gbs', 'GHF')  # -147.634028141 a.u.
    # start('CO.mol', 'STO-3G.gbs', 'MP2')  # -111.354512528 a.u.
    # start('H2O.mol', 'STO-3G.gbs', 'RHF', True)
    start('C2H4.mol', '3-21G.gbs', 'RHF', True)  # -77.600460844 a.u 20.19939799526804s


def start(mol, basis, method, symmetry=False):
    np.set_printoptions(linewidth=100000, threshold=np.inf)
    start_time = time.clock()
    electron_energy = correlation = 0.0

    nuclei_array, electrons, multiplicity = read_mol_file(mol)

    if symmetry:
        molecule = MoleculeFactory().point_group(nuclei_array)
    else:
        molecule = Molecule(nuclei_array, None)

    basis_set_array = read_basis_set_file(basis, molecule.nuclei_array)
    symmetry_object = Symmetry(molecule, basis_set_array)
    print(symmetry_object.symmetry_matrix if symmetry else '', end='\n\n')

    coulomb_law_matrix = coulomb_matrix(nuclei_array)
    nuclear_repulsion = coulomb_law_matrix.sum() / 2

    print('*******************************************************************************************************')
    print('\nA BASIC QUANTUM CHEMICAL PROGRAM IN PYTHON\n\n')
    print([x.element for x in nuclei_array])
    print('\nNUCLEAR REPULSION ARRAY')
    print(coulomb_law_matrix)

    if method == 'RHF':
        electron_energy = RestrictedHF(molecule.nuclei_array, basis_set_array, electrons, symmetry_object).begin()[0]
    if method == 'UHF':
        electron_energy = DODSUnrestricted(molecule.nuclei_array, basis_set_array, electrons, multiplicity,
        symmetry_object).begin()[0]
    if method == 'CUHF':
        electron_energy = ConstrainedUnrestricted(molecule.nuclei_array, basis_set_array, electrons, multiplicity,
        symmetry_object).begin()[0]
    if method == 'GHF':
        electron_energy = BlockedHartreeFock(molecule.nuclei_array, basis_set_array, electrons, multiplicity,
        symmetry_object).begin()[0]
    if method == 'MP2':
        electron_energy, correlation = MoellerPlesset.second_order(molecule.nuclei_array, basis_set_array, electrons,
        symmetry_object)

    print('NUCLEAR REPULSION ENERGY:    ' + str(nuclear_repulsion) + ' a.u.')
    print('SCF ENERGY:                  ' + str(electron_energy) + ' a.u.')
    print('CORRELATION ENERGY:          ' + str(correlation) + ' a.u.')
    print('TOTAL ENERGY:                ' + str(electron_energy + nuclear_repulsion + correlation) + ' a.u.')
    print('\n*****************************************************************************************************')
    print('\nTIME TAKEN: ' + str(time.clock() - start_time) + 's')
    print("\nWhat I cannot create I cannot understand - Richard Feynman\n")


if __name__ == "__main__":
    menu()
