import numpy as np
import time
from src.main.matrixelements import GMatrixElementUnrestricted as GMatrixElement
from src.main.matrixelements import DensityMatrixUnrestricted as DensityMatrixElement
from src.main.common import TotalEnergy
from math import floor, ceil


class UHFSCFProcedure:

    def __init__(self, core_hamiltonian_matrix, transformation_matrix, matrix, electrons, repulsion_dictionary):
        self.core_hamiltonian_matrix = core_hamiltonian_matrix
        self.transformation_matrix = transformation_matrix
        self.matrix = matrix
        self.electrons = electrons
        self.repulsion = repulsion_dictionary

    def begin_scf(self, orbital_coefficients, multiplicity):
        start = time.clock()
        difference = floor(multiplicity / 2)
        electron_alpha = ceil(self.electrons / 2) + difference
        electron_beta = floor(self.electrons / 2) - difference
        orbital_coefficients_alpha = orbital_coefficients
        orbital_coefficients_beta = orbital_coefficients
        total_energy = 0
        previous_total_energy = 0
        delta_energy = 1
        eigenvalues_alpha = []
        eigenvalues_beta = []

        while abs(delta_energy) > 1e-9:

            if total_energy == 0:
                density_matrix_alpha = self.matrix.create_matrix(DensityMatrixElement(orbital_coefficients_alpha, (electron_alpha + 1)))
                density_matrix_beta = self.matrix.create_matrix(DensityMatrixElement(orbital_coefficients_beta, (electron_beta - 1)))
            else:
                density_matrix_alpha = self.matrix.create_matrix(DensityMatrixElement(orbital_coefficients_alpha, electron_alpha))
                density_matrix_beta = self.matrix.create_matrix(DensityMatrixElement(orbital_coefficients_beta, electron_beta))
            density_matrix_total = density_matrix_alpha + density_matrix_beta

            g_matrix_alpha = self.matrix.create_matrix(GMatrixElement(density_matrix_total, density_matrix_alpha, self.repulsion))
            g_matrix_beta = self.matrix.create_matrix(GMatrixElement(density_matrix_total, density_matrix_beta, self.repulsion))
            fock_matrix_alpha = self.core_hamiltonian_matrix + g_matrix_alpha
            fock_matrix_beta = self.core_hamiltonian_matrix + g_matrix_beta

            total_energy = TotalEnergy.unrestricted(density_matrix_total, density_matrix_alpha, density_matrix_beta, self.core_hamiltonian_matrix, fock_matrix_alpha, fock_matrix_beta)
            delta_energy = previous_total_energy - total_energy
            previous_total_energy = total_energy
            print('SCF ENERGY: ' + str(total_energy) + ' a.u.')

            orthonormal_fock_matrix = np.transpose(self.transformation_matrix) * fock_matrix_alpha * self.transformation_matrix
            eigenvalues, eigenvectors = np.linalg.eigh(orthonormal_fock_matrix)
            sort = np.argsort(eigenvalues)
            eigenvalues_alpha = np.array(eigenvalues)[sort]
            eigenvectors = eigenvectors[:, sort]
            orbital_coefficients_alpha = self.transformation_matrix * eigenvectors

            orthonormal_fock_matrix = np.transpose(self.transformation_matrix) * fock_matrix_beta * self.transformation_matrix
            eigenvalues, eigenvectors = np.linalg.eigh(orthonormal_fock_matrix)
            sort = np.argsort(eigenvalues)
            eigenvalues_beta = np.array(eigenvalues)[sort]
            eigenvectors = eigenvectors[:, sort]
            orbital_coefficients_beta = self.transformation_matrix * eigenvectors

        print('TIME TAKEN: ' + str(time.clock() - start) + 's\n')

        print('\nORBITAL ENERGY EIGENVALUES')
        print(eigenvalues_alpha)
        print(eigenvalues_beta)

        print('\nORBITAL COEFFICIENTS')
        print(orbital_coefficients_alpha, end='\n\n')
        print(orbital_coefficients_beta, end='\n\n')

        return total_energy
