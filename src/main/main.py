from src.main.hartreefock import UnrestrictedHF, RestrictedHF
from src.main.fileinput import FileInputBasis, FileInputNuclei
from src.main.common import CoulombsLawMatrix
from src.main.moellerplesset import MoellerPlesset
import numpy as np
import time


def menu():
    # start('CO.mol', 'STO-3G.gbs', 'MP2')
    # start('O2.mol', 'STO-3G.gbs', 'UHF')
    start('C2H4.mol', '3-21G.gbs', 'MP2')


def start(mol, basis, method):
    np.set_printoptions(linewidth=100000, threshold=np.inf)
    start_time = time.clock()
    electron_energy = correlation = 0.0

    nuclei_array, electrons, multiplicity = FileInputNuclei.read(mol)
    basis_set_array = FileInputBasis.read(basis, nuclei_array)

    coulomb_law_matrix = CoulombsLawMatrix.create(nuclei_array)
    nuclear_repulsion = coulomb_law_matrix.sum() / 2

    print('*******************************************************************************************************')
    print('\nA BASIC QUANTUM CHEMICAL PROGRAM IN PYTHON\n\n')
    print([x.element for x in nuclei_array])
    print('\nNUCLEAR REPULSION ARRAY')
    print(coulomb_law_matrix)

    if method == 'RHF':
        electron_energy = RestrictedHF(nuclei_array, basis_set_array, electrons).begin()[0]
    if method == 'UHF':
        electron_energy = UnrestrictedHF(nuclei_array, basis_set_array, electrons, multiplicity).begin()[0]
    if method == 'MP2':
        electron_energy, correlation = MoellerPlesset.second_order(nuclei_array, basis_set_array, electrons)

    print('NUCLEAR REPULSION ENERGY:    ' + str(nuclear_repulsion) + ' a.u.')
    print('SCF ENERGY:                  ' + str(electron_energy) + ' a.u.')
    print('CORRELATION ENERGY:          ' + str(correlation) + ' a.u.')
    print('TOTAL ENERGY:                ' + str(electron_energy + nuclear_repulsion + correlation) + ' a.u.')
    print('\n*****************************************************************************************************')
    print('\nTIME TAKEN: ' + str(time.clock() - start_time) + 's')
    print("\nWhat I cannot create I cannot understand - Richard Feynman\n")


if __name__ == "__main__":
    menu()
