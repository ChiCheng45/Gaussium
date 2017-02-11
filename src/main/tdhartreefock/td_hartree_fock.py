from src.main.matrixelements import molecular_orbitals
from src.main.matrixelements import spin_basis_set
from src.main.matrixelements import spin_orbital_energies
from src.main.tdhartreefock import TDHFMatrixSymmetryRestricted
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

        tdhf = TDHFMatrixSymmetryRestricted(repulsion, orbital_energies, occupied_orbitals, unoccupied_orbitals)
        singlet_a_matrix, singlet_b_matrix = tdhf.create_matrices_singlet()
        triplet_a_matrix, triplet_b_matrix = tdhf.create_matrices_triplet()

        if not tda:
            print('BEGIN TDHF EXCITED STATE CALCULATION\n')
            singlet_tdhf_matrix = (singlet_a_matrix + singlet_b_matrix) * (singlet_a_matrix - singlet_b_matrix)
            singlet_excitation_energies = np.linalg.eig(singlet_tdhf_matrix)[0]
            singlet_excitation_energies = [sqrt(x) for x in singlet_excitation_energies]
            singlet_excitation_energies.sort()
            triplet_tdhf_matrix = (triplet_a_matrix + triplet_b_matrix) * (triplet_a_matrix - triplet_b_matrix)
            triplet_excitation_energies = np.linalg.eig(triplet_tdhf_matrix)[0]
            triplet_excitation_energies = [sqrt(x) for x in triplet_excitation_energies]
            triplet_excitation_energies.sort()
        else:
            print('BEGIN CIS EXCITED STATE CALCULATION\n')
            singlet_excitation_energies = np.linalg.eigh(singlet_a_matrix)[0]
            triplet_excitation_energies = np.linalg.eigh(triplet_a_matrix)[0]

        print('SINGLET EXCITATION ENERGIES\n{}\n'.format(singlet_excitation_energies))
        print('TRIPLET EXCITATION ENERGIES\n{}\n\n'.format(triplet_excitation_energies))

        return electron_energy, 0.0
