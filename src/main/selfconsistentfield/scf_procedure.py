import numpy as np

from src.main.matrixelements.g_matrix import GMatrixElements
from src.main.matrixelements.density_matrix import DensityMatrixElement


class SCFProcedure:

    def __init__(self, core_hamiltonian_matrix, transformation_matrix, matrix, total_energy_calc, basis_set_array, electrons, orbital_coefficients, orbital_energy_matrix):
        self.core_hamiltonian_matrix = core_hamiltonian_matrix
        self.transformation_matrix = transformation_matrix
        self.matrix = matrix
        self.total_energy_calc = total_energy_calc
        self.basis_set_array = basis_set_array
        self.orbital_coefficients = orbital_coefficients
        self.electrons = electrons
        self.orbital_energy_matrix = orbital_energy_matrix

    def begin_scf(self):

        total_energy = 0
        previous_total_energy = 0
        iteration_counter = 0
        delta_energy = 1

        while abs(delta_energy) > 0.00000000001:

            iteration_counter += 1
            print('\n\n----------ITERATION: ' + str(iteration_counter) + str('----------'))

            print('\nORBITAL ENERGY EIGENVALUES')
            print(self.orbital_energy_matrix)

            print('\nORBITAL COEFFICIENTS')
            print(self.orbital_coefficients)

            density_matrix_element = DensityMatrixElement(self.orbital_coefficients, self.electrons)
            density_matrix = self.matrix.create_matrix(density_matrix_element)
            print('\nDENSITY MATRIX')
            print(density_matrix)

            g_matrix_elements = GMatrixElements(density_matrix, self.basis_set_array)
            g_matrix = self.matrix.create_matrix(g_matrix_elements)
            fock_matrix = self.core_hamiltonian_matrix + g_matrix
            print('\nFOCK MATRIX')
            print(fock_matrix)

            total_energy = self.total_energy_calc.calculate_total_energy(density_matrix, self.core_hamiltonian_matrix, fock_matrix)
            delta_energy = previous_total_energy - total_energy
            previous_total_energy = total_energy
            print('\nSCF ENERGY: ' + str(total_energy) + ' a.u.')

            orthonormal_fock_matrix = self.transformation_matrix.T * fock_matrix * self.transformation_matrix
            eigenvalues, eigenvectors = np.linalg.eig(orthonormal_fock_matrix)
            sort = eigenvalues.argsort()[::1]
            eigenvalues = eigenvalues[sort]
            eigenvectors = eigenvectors[:, sort]

            self.orbital_energy_matrix = np.diag(eigenvalues)
            self.orbital_coefficients = self.transformation_matrix * eigenvectors

        return total_energy
