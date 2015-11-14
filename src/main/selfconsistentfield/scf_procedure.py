import numpy as np
from src.main.matrixelements.g_matrix_element import GMatrixElement
from src.main.matrixelements.density_matrix_element import DensityMatrixElement
from src.main.selfconsistentfield import TotalEnergy


class SCFProcedure:

    def __init__(self, core_hamiltonian_matrix, transformation_matrix, matrix, basis_set_array, electrons):
        self.core_hamiltonian_matrix = core_hamiltonian_matrix
        self.transformation_matrix = transformation_matrix
        self.matrix = matrix
        self.basis_set_array = basis_set_array
        self.electrons = electrons

    def begin_scf(self, orbital_coefficients, repulsion_dictionary):

        total_energy = 0
        previous_total_energy = 0
        delta_energy = 1

        while abs(delta_energy) > 1e-9:

            density_matrix_element = DensityMatrixElement(orbital_coefficients, self.electrons)
            density_matrix = self.matrix.create_matrix(density_matrix_element)

            g_matrix_elements = GMatrixElement(density_matrix, repulsion_dictionary)
            g_matrix = self.matrix.create_matrix(g_matrix_elements)
            fock_matrix = self.core_hamiltonian_matrix + g_matrix

            total_energy = TotalEnergy.calculate_total_energy(density_matrix, self.core_hamiltonian_matrix, fock_matrix)
            delta_energy = previous_total_energy - total_energy
            previous_total_energy = total_energy
            print('SCF ENERGY: ' + str(total_energy) + ' a.u.')

            orthonormal_fock_matrix = self.transformation_matrix.T * fock_matrix * self.transformation_matrix
            eigenvalues, eigenvectors = np.linalg.eig(orthonormal_fock_matrix)
            sort = eigenvalues.argsort()[::1]
            eigenvalues = eigenvalues[sort]
            eigenvectors = eigenvectors[:, sort]

            orbital_energy_matrix = np.diag(eigenvalues)
            orbital_coefficients = self.transformation_matrix * eigenvectors

        print('\nORBITAL ENERGY EIGENVALUES')
        print(orbital_energy_matrix)

        print('\nORBITAL COEFFICIENTS')
        print(orbital_coefficients, end='\n\n')

        return total_energy
