from src.main.matrixelements import molecular_orbitals
from src.main.matrixelements import spin_basis_set
from src.main.matrixelements import spin_orbital_energies
from src.main.tdhartreefock import TDHFMatrix
from math import sqrt
import numpy as np


class TimeDependentHartreeFock:

    def __init__(self, hartree_fock):
        self.hartree_fock = hartree_fock

    def calculate(self, tda=False):
        electron_energy, orbital_energies, orbital_coefficients = self.hartree_fock.begin_scf()

        repulsion = spin_basis_set(molecular_orbitals(self.hartree_fock.repulsion, orbital_coefficients))
        orbital_energies = spin_orbital_energies(orbital_energies)
        occupied_orbitals = self.hartree_fock.electrons
        unoccupied_orbitals = len(orbital_energies) - occupied_orbitals

        a_matrix, b_matrix = TDHFMatrix(repulsion, orbital_energies, occupied_orbitals,
        unoccupied_orbitals).create_matrices()

        if not tda:
            print('BEGIN TDHF EXCITED STATE CALCULATION\n')
            tdhf_matrix = (a_matrix + b_matrix) * (a_matrix - b_matrix)
            excitation_energies = np.linalg.eig(tdhf_matrix)[0]
            excitation_energies = [sqrt(x) for x in excitation_energies]
            excitation_energies.sort()
        else:
            print('BEGIN CIS EXCITED STATE CALCULATION\n')
            excitation_energies = np.linalg.eigh(a_matrix)[0]

        print('EXCITATION ENERGIES\n{}\n\n'.format(excitation_energies))

        return electron_energy, orbital_energies, orbital_coefficients
