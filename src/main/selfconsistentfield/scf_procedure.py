import numpy as np

from src.main.matrixelements.g_matrix import GMatrixElements
from src.main.matrixelements.density_matrix import DensityMatrixElement


class SCFProcedure:

    """
    To begin the scf procedure we will need the initial guess for the density matrix, the two electron repulsion
    matrixelements, the kinetic and nuclear potential energy matrix to form the H_core, the transformation matrix and the
    basis_set_array.
    """
    def __init__(self, core_hamiltonian_matrix, transformation_matrix, matrix, total_energy, nuclear_repulsion, basis_set_array, electrons):
        self.core_hamiltonian_matrix = core_hamiltonian_matrix
        self.transformation_matrix = transformation_matrix
        self.matrix = matrix
        self.total_energy = total_energy
        self.nuclear_repulsion = nuclear_repulsion
        self.basis_set_array = basis_set_array
        self.electrons = electrons
        self.previous_total_energy = 0
        self.delta_energy = 0
        self.iteration_counter = 0

    def begin_scf(self, density_matrix):

        self.iteration_counter += 1
        print('\n\n----------ITERATION: ' + str(self.iteration_counter) + str('----------'))

        print('\nDENSITY MATRIX')
        print(density_matrix)

        g_matrix_elements = GMatrixElements(density_matrix, self.basis_set_array)
        g_matrix = self.matrix.create_matrix(g_matrix_elements)
        fock_matrix = self.core_hamiltonian_matrix + g_matrix
        print('\nFOCK MATRIX')
        print(fock_matrix)

        orthonormal_fock_matrix = self.transformation_matrix.T * fock_matrix * self.transformation_matrix
        eigenvalues,eigenvectors = np.linalg.eig(orthonormal_fock_matrix)
        sort = eigenvalues.argsort()[::1]
        eigenvalues = eigenvalues[sort]
        eigenvectors = eigenvectors[:, sort]

        orbital_energy_matrix = np.diag(eigenvalues)
        print('\nORBITAL ENERGY EIGENVALUES')
        print(orbital_energy_matrix)

        orbital_coefficients = self.transformation_matrix * eigenvectors
        print('\nORBITAL COEFFICIENTS')
        print(orbital_coefficients)

        e_total = self.total_energy.calculate_total_energy(density_matrix, self.core_hamiltonian_matrix, fock_matrix)
        print('\nELECTRON ENERGY: ' + str(e_total) + ' a.u.')
        print('TOTAL ENERGY: ' + str(e_total + self.nuclear_repulsion) + ' a.u.')

        self.delta_energy = self.previous_total_energy - e_total
        self.previous_total_energy = e_total

        """
        Form the density matrix for the next iteration and solve the for the total_energy recursively this is easier
        to implement but it might cause a stack overflow when we get into bigger basis set and bigger molecules. This is
        also probably inefficient.
        """
        densitymatrix = DensityMatrixElement(orbital_coefficients, self.electrons)
        density_matrix = self.matrix.create_matrix(densitymatrix)

        if abs(self.delta_energy) > 0.00000000001:
            return self.begin_scf(density_matrix)
        else:
            return e_total
