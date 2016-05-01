import numpy as np


class LinearAlgebra:

    def __init__(self, orbital_overlap):
        self.transformation_matrix = self.create_transformation_matrix(orbital_overlap)

    def diagonalize(self, fock_matrix):
        orthonormal_h_matrix = self.transformation_matrix.T * fock_matrix * self.transformation_matrix
        orbital_energies, orbital_coefficients = np.linalg.eigh(orthonormal_h_matrix)
        sort = np.argsort(orbital_energies)
        orbital_energies = np.array(orbital_energies)[sort]
        orbital_coefficients = orbital_coefficients[:, sort]
        orbital_coefficients = self.transformation_matrix * orbital_coefficients
        return orbital_energies, orbital_coefficients

    def create_transformation_matrix(self, orbital_overlap):
        s_matrix_eigenvalues, s_matrix_unitary = np.linalg.eigh(orbital_overlap)
        sort = np.argsort(s_matrix_eigenvalues)
        s_matrix_eigenvalues = np.array(s_matrix_eigenvalues)[sort]
        s_matrix_unitary = s_matrix_unitary[:, sort]
        s_matrix_eigenvalues = [x**(-1/2) for x in s_matrix_eigenvalues]
        x_canonical = s_matrix_unitary * np.diag(s_matrix_eigenvalues)
        return x_canonical

    def orthonormalize(self, matrix):
        return self.transformation_matrix.T * matrix * self.transformation_matrix

    def non_orthonormal(self, matrix):
        return np.linalg.inv(self.transformation_matrix.T) * matrix * np.linalg.inv(self.transformation_matrix)


class BlockedLinearAlgebra(LinearAlgebra):

    def diagonalize(self, fock_matrix):
        orthonormal_h_matrix = self.transformation_matrix.T * fock_matrix * self.transformation_matrix
        orbital_energies, orbital_coefficients = np.linalg.eigh(orthonormal_h_matrix)

        sort = np.argsort(np.fabs(orbital_coefficients.tolist()[0]))[::-1]
        orbital_energies = np.array(orbital_energies)[sort]
        orbital_coefficients = orbital_coefficients[:, sort]

        half_matrix_size = orbital_energies.shape[0] // 2
        orbital_energies_alph = orbital_energies[:half_matrix_size]
        orbital_energies_beta = orbital_energies[half_matrix_size:]

        sort_alph = np.argsort(orbital_energies_alph)
        sort_beta = [i + half_matrix_size for i in np.argsort(orbital_energies_beta)]
        sort = np.concatenate((sort_alph, sort_beta))

        orbital_energies = np.array(orbital_energies)[sort]
        orbital_coefficients = orbital_coefficients[:, sort]

        orbital_coefficients = self.transformation_matrix * orbital_coefficients
        return orbital_energies, orbital_coefficients

    def create_transformation_matrix(self, orbital_overlap):
        s_matrix_eigenvalues, s_matrix_unitary = np.linalg.eigh(orbital_overlap)

        sort = np.argsort(np.fabs(s_matrix_unitary.tolist()[0]))[::-1]
        s_matrix_eigenvalues = np.array(s_matrix_eigenvalues)[sort]
        s_matrix_unitary = s_matrix_unitary[:, sort]

        half_matrix_size = s_matrix_eigenvalues.shape[0] // 2
        s_matrix_eigenvalues_alph = s_matrix_eigenvalues[:half_matrix_size]
        s_matrix_eigenvalues_beta = s_matrix_eigenvalues[half_matrix_size:]

        sort_alph = np.argsort(s_matrix_eigenvalues_alph)
        sort_beta = [i + half_matrix_size for i in np.argsort(s_matrix_eigenvalues_beta)]
        sort = np.concatenate((sort_alph, sort_beta))

        s_matrix_eigenvalues = np.array(s_matrix_eigenvalues)[sort]
        s_matrix_unitary = s_matrix_unitary[:, sort]
        s_matrix_eigenvalues = [x**(-1/2) for x in s_matrix_eigenvalues]
        x_canonical = s_matrix_unitary * np.diag(s_matrix_eigenvalues)
        return x_canonical
