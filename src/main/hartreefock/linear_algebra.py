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
        half_matrix_size = fock_matrix.shape[0] // 2
        zeros = np.zeros((half_matrix_size, half_matrix_size))

        top_left = orthonormal_h_matrix[:half_matrix_size, :half_matrix_size]
        bot_righ = orthonormal_h_matrix[half_matrix_size:, half_matrix_size:]
        orbital_energies_alph, orbital_coefficients_alph = np.linalg.eigh(top_left)
        orbital_energies_beta, orbital_coefficients_beta = np.linalg.eigh(bot_righ)

        sort = np.argsort(orbital_energies_alph)
        orbital_energies_alph = np.array(orbital_energies_alph)[sort]
        orbital_coefficients_alph = orbital_coefficients_alph[:, sort]

        sort = np.argsort(orbital_energies_beta)
        orbital_energies_beta = np.array(orbital_energies_beta)[sort]
        orbital_coefficients_beta = orbital_coefficients_beta[:, sort]


        orbital_coefficients = np.bmat([
                [orbital_coefficients_alph, zeros],
                [zeros, orbital_coefficients_beta]
        ])

        orbital_energies = np.append(orbital_energies_alph, orbital_energies_beta)
        orbital_coefficients = self.transformation_matrix * orbital_coefficients
        return orbital_energies, orbital_coefficients

    def create_transformation_matrix(self, orbital_overlap):
        half_matrix_size = orbital_overlap.shape[0] // 2
        zeros = np.zeros((half_matrix_size, half_matrix_size))

        top_left = orbital_overlap[:half_matrix_size, :half_matrix_size]
        bot_righ = orbital_overlap[half_matrix_size:, half_matrix_size:]
        s_matrix_eigenvalues_top, s_matrix_unitary_top = np.linalg.eigh(top_left)
        s_matrix_eigenvalues_bot, s_matrix_unitary_bot = np.linalg.eigh(bot_righ)

        sort = np.argsort(s_matrix_eigenvalues_top)
        s_matrix_eigenvalues_top = np.array(s_matrix_eigenvalues_top)[sort]
        s_matrix_unitary_top = s_matrix_unitary_top[:, sort]

        sort = np.argsort(s_matrix_eigenvalues_bot)
        s_matrix_eigenvalues_bot = np.array(s_matrix_eigenvalues_bot)[sort]
        s_matrix_unitary_bot = s_matrix_unitary_bot[:, sort]

        s_matrix_unitary = np.bmat([
                [s_matrix_unitary_top, zeros],
                [zeros, s_matrix_unitary_bot]
        ])

        s_matrix_eigenvalues = np.append(s_matrix_eigenvalues_top, s_matrix_eigenvalues_bot)
        s_matrix_eigenvalues = [x**(-1/2) for x in s_matrix_eigenvalues]
        x_canonical = s_matrix_unitary * np.diag(s_matrix_eigenvalues)
        return x_canonical
