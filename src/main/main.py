from src.main.hartreefock import HartreeFock
from src.main.fileinput import FileInputBasis, FileInputNuclei
from src.main.common import CoulombsLawArray
from src.main.moellerplesset import MoellerPlesset
import numpy as np
import time


def menu():
    start('C2H4.mol', '3-21G.gbs', 'RHF')
    start('O2.mol', 'STO-3G.gbs', 'UHF')
    start('HeH+.mol', 'STO-3g.gbs', 'MP2')


def start(mol, basis, method):
    np.set_printoptions(linewidth=100000, threshold=np.inf)
    start_time = time.clock()
    electron_energy = correlation = 0.0

    nuclei_array, electrons, multiplicity = FileInputNuclei.read(mol)
    basis_set_array = FileInputBasis.read(basis, nuclei_array)
    nuclei_name_list = [x.element for x in nuclei_array]

    coulomb_law_array = CoulombsLawArray.calculate_total_electric_potential_energy(nuclei_array)
    nuclear_repulsion = coulomb_law_array.sum() / 2

    print('*******************************************************************************************************')
    print('\nA BASIC QUANTUM CHEMICAL PROGRAM IN PYTHON\n\n')
    print(nuclei_name_list)
    print('\nNUCLEAR REPULSION ARRAY')
    print(coulomb_law_array)

    if method == 'RHF':
        electron_energy = HartreeFock(nuclei_array, electrons, multiplicity, basis_set_array).restricted()[0]
    if method == 'UHF':
        electron_energy = HartreeFock(nuclei_array, electrons, multiplicity, basis_set_array).unrestricted()[0]
    if method == 'MP2':
        electron_energy, correlation = MoellerPlesset.second_order(nuclei_array, electrons, multiplicity, basis_set_array)

    print('\nNUCLEAR REPULSION ENERGY:    ' + str(nuclear_repulsion) + ' a.u.')
    print('SCF ENERGY:                  ' + str(electron_energy) + ' a.u.')
    print('CORRELATION ENERGY:          ' + str(correlation) + ' a.u.')
    print('TOTAL ENERGY:                ' + str(electron_energy + nuclear_repulsion + correlation) + ' a.u.')
    print('\n*****************************************************************************************************')
    print('\nTIME TAKEN: ' + str(time.clock() - start_time) + 's')
    print("\nWhat I cannot create I cannot understand - Richard Feynman\n")


if __name__ == "__main__":
    menu()
