import numpy as np
from src.main.matrixelements.g_matrix_element import GMatrixElement
from src.main.matrixelements.density_matrix_element import DensityMatrixElement
from src.main.selfconsistentfield import TotalEnergy


class SCFProcedure:

    def __init__(self, core_hamiltonian_matrix, transformation_matrix, matrix, electrons, repulsion_dictionary):
        self.core_hamiltonian_matrix = core_hamiltonian_matrix
        self.transformation_matrix = transformation_matrix
        self.matrix = matrix
        self.electrons = electrons
        self.repulsion_dictionary = repulsion_dictionary

    def begin_scf(self, orbital_coefficients):
        total_energy = 0
        previous_total_energy = 0
        delta_energy = 1
        eigenvalues = []

        while abs(delta_energy) > 1e-9:
            density_matrix = self.matrix.create_matrix(DensityMatrixElement(orbital_coefficients, self.electrons))
            g_matrix = self.matrix.create_matrix(GMatrixElement(density_matrix, self.repulsion_dictionary))
            fock_matrix = self.core_hamiltonian_matrix + g_matrix

            total_energy = TotalEnergy.calculate_total_energy(density_matrix, self.core_hamiltonian_matrix, fock_matrix)
            delta_energy = previous_total_energy - total_energy
            previous_total_energy = total_energy
            print('SCF ENERGY: ' + str(total_energy) + ' a.u.')

            orthonormal_fock_matrix = np.transpose(self.transformation_matrix) * fock_matrix * self.transformation_matrix
            eigenvalues, eigenvectors = np.linalg.eigh(orthonormal_fock_matrix)
            sort = np.argsort(eigenvalues)
            eigenvalues = np.array(eigenvalues)[sort]
            eigenvectors = eigenvectors[:, sort]

            orbital_coefficients = self.transformation_matrix * eigenvectors

        print('\nORBITAL ENERGY EIGENVALUES')
        print(eigenvalues)

        print('\nORBITAL COEFFICIENTS')
        print(orbital_coefficients, end='\n\n')

        return total_energy
