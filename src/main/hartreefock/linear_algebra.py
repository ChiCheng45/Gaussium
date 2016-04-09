import numpy as np


class LinearAlgebra:

    def __init__(self, orbital_overlap):
        self.transformation_matrix = self.create_transformation_matrix(orbital_overlap)

    def diagonalize(self, fock_matrix):
        orthonormal_h_matrix = np.transpose(self.transformation_matrix) * fock_matrix * self.transformation_matrix
        orbital_energies, orbital_coefficients = np.linalg.eigh(orthonormal_h_matrix)
        sort = np.argsort(orbital_energies)
        orbital_energies = np.array(orbital_energies)[sort]
        orbital_coefficients = orbital_coefficients[:, sort]
        orbital_coefficients = self.transformation_matrix * orbital_coefficients
        return orbital_energies, orbital_coefficients

    def orthonormalize(self, matrix):
        return np.transpose(self.transformation_matrix) * matrix * self.transformation_matrix

    def non_orthonormal(self, matrix):
        return np.linalg.inv(np.transpose(self.transformation_matrix)) * matrix * np.linalg.inv(self.transformation_matrix)

    @staticmethod
    def create_transformation_matrix(orbital_overlap):
        s_matrix_eigenvalues, s_matrix_unitary = np.linalg.eigh(orbital_overlap)
        sort = np.argsort(s_matrix_eigenvalues)
        s_matrix_eigenvalues = np.array(s_matrix_eigenvalues)[sort]
        s_matrix_unitary = s_matrix_unitary[:, sort]
        s_matrix_eigenvalues = [x**(-1/2) for x in s_matrix_eigenvalues]
        x_canonical = s_matrix_unitary * np.diag(s_matrix_eigenvalues)
        return x_canonical
