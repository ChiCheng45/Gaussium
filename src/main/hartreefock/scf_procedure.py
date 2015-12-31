from src.main.matrixelements import GElementUnrestricted, GElementRestricted
from src.main.matrixelements import DensityElementUnrestricted, DensityElementRestricted
from src.main.common import TotalEnergy
import numpy as np
from math import floor, ceil
import time


class SCF:

    def __init__(self, core_hamiltonian, transformation_matrix, create_matrix, electrons, repulsion_dictionary):
        self.core_hamiltonian = core_hamiltonian
        self.transformation_matrix = transformation_matrix
        self.create = create_matrix
        self.electrons = electrons
        self.repulsion = repulsion_dictionary
        self.total_energy = 0
        self.previous_total_energy = 0
        self.delta_energy = 1

    def begin_rhf(self, orbital_coefficients):
        start = time.clock()
        eigenvalues = []

        while abs(self.delta_energy) > 1e-9:
            density_matrix = self.create(DensityElementRestricted(orbital_coefficients, self.electrons))
            g_matrix = self.create(GElementRestricted(density_matrix, self.repulsion))
            fock_matrix = self.core_hamiltonian + g_matrix
            eigenvalues, orbital_coefficients = self.diagonalize(fock_matrix)

            self.total_energy = TotalEnergy.restricted(density_matrix, self.core_hamiltonian, fock_matrix)
            self.delta_energy = self.previous_total_energy - self.total_energy
            self.previous_total_energy = self.total_energy
            print('SCF ENERGY: ' + str(self.total_energy) + ' a.u.')

        print('TIME TAKEN: ' + str(time.clock() - start) + 's\n')
        print('\nORBITAL ENERGY EIGENVALUES')
        print(eigenvalues)
        print('\nORBITAL COEFFICIENTS')
        print(orbital_coefficients, end='\n\n')
        return self.total_energy

    def begin_uhf(self, orbital_coefficients, multiplicity):
        start = time.clock()
        difference = floor(multiplicity / 2)
        electrons_alpha = ceil(self.electrons / 2) + difference
        electrons_beta = floor(self.electrons / 2) - difference
        orbital_coefficients_alpha = orbital_coefficients
        orbital_coefficients_beta = orbital_coefficients
        eigenvalues_alpha = []
        eigenvalues_beta = []

        while abs(self.delta_energy) > 1e-9:
            if self.total_energy == 0:
                density_matrix_alpha = self.create(DensityElementUnrestricted(orbital_coefficients_alpha, (electrons_alpha + 1)))
                density_matrix_beta = self.create(DensityElementUnrestricted(orbital_coefficients_beta, (electrons_beta - 1)))
            else:
                density_matrix_alpha = self.create(DensityElementUnrestricted(orbital_coefficients_alpha, electrons_alpha))
                density_matrix_beta = self.create(DensityElementUnrestricted(orbital_coefficients_beta, electrons_beta))

            density_matrix_total = density_matrix_alpha + density_matrix_beta
            g_matrix_alpha = self.create(GElementUnrestricted(density_matrix_total, density_matrix_alpha, self.repulsion))
            g_matrix_beta = self.create(GElementUnrestricted(density_matrix_total, density_matrix_beta, self.repulsion))
            fock_matrix_alpha = self.core_hamiltonian + g_matrix_alpha
            fock_matrix_beta = self.core_hamiltonian + g_matrix_beta
            eigenvalues_alpha, orbital_coefficients_alpha = self.diagonalize(fock_matrix_alpha)
            eigenvalues_beta, orbital_coefficients_beta = self.diagonalize(fock_matrix_beta)

            self.total_energy = TotalEnergy.unrestricted(density_matrix_total, density_matrix_alpha, density_matrix_beta, self.core_hamiltonian, fock_matrix_alpha, fock_matrix_beta)
            self.delta_energy = self.previous_total_energy - self.total_energy
            self.previous_total_energy = self.total_energy
            print('SCF ENERGY: ' + str(self.total_energy) + ' a.u.')

        print('TIME TAKEN: ' + str(time.clock() - start) + 's\n')
        print('\nALPHA ORBITAL ENERGY EIGENVALUES')
        print(eigenvalues_alpha)
        print('\nBETA ORBITAL ENERGY EIGENVALUES')
        print(eigenvalues_beta)
        print('\nALPHA ORBITAL COEFFICIENTS')
        print(orbital_coefficients_alpha, end='\n\n')
        print('\nBETA ORBITAL COEFFICIENTS')
        print(orbital_coefficients_beta, end='\n\n')
        return self.total_energy

    def diagonalize(self, fock_matrix):
        orthonormal_fock_matrix = np.transpose(self.transformation_matrix) * fock_matrix * self.transformation_matrix
        eigenvalues, eigenvectors = np.linalg.eigh(orthonormal_fock_matrix)
        sort = np.argsort(eigenvalues)
        eigenvalues = np.array(eigenvalues)[sort]
        eigenvectors = eigenvectors[:, sort]
        orbital_coefficients = self.transformation_matrix * eigenvectors
        return eigenvalues, orbital_coefficients
