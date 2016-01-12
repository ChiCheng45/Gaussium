from src.main.matrixelements import Matrix
import numpy as np


class GMatrixRestricted(Matrix):

    def __init__(self, repulsion_dictionary):
        super().__init__()
        self.repulsion_dictionary = repulsion_dictionary
        self.density_matrix = np.matrix([])

    def create(self, density_matrix):
        self.density_matrix = density_matrix
        self.matrix_size = density_matrix.shape[0]
        return self.create_matrix(self.calculate_restricted)

    def calculate_restricted(self, i, j):
        g_ij = 0
        for a in range(self.matrix_size):
            for b in range(self.matrix_size):
                coulomb_integral = self.repulsion_dictionary.item(i, j, a, b)
                exchange_integral = self.repulsion_dictionary.item(i, b, a, j)
                g_ij += self.density_matrix.item(a, b) * (coulomb_integral - (1/2) * exchange_integral)
        return g_ij


class GMatrixUnrestricted(Matrix):

    def __init__(self, repulsion_dictionary):
        super().__init__()
        self.repulsion_dictionary = repulsion_dictionary
        self.density_matrix_total = np.matrix([])
        self.density_matrix_alph = np.matrix([])
        self.density_matrix_beta = np.matrix([])
        self.density_matrix = np.matrix([])

    def set_density_matrices(self, density_matrix_alph, density_matrix_beta):
        self.matrix_size = density_matrix_alph.shape[0]
        self.density_matrix_alph = density_matrix_alph
        self.density_matrix_beta = density_matrix_beta
        self.density_matrix_total = density_matrix_alph + density_matrix_beta

    def create_alph(self):
        self.density_matrix = self.density_matrix_alph
        return self.create_matrix(self.calculate_unrestricted)

    def create_beta(self):
        self.density_matrix = self.density_matrix_beta
        return self.create_matrix(self.calculate_unrestricted)

    def calculate_unrestricted(self, i, j):
        g_ij = 0
        for a in range(self.matrix_size):
            for b in range(self.matrix_size):
                coulomb_integral = self.repulsion_dictionary.item(i, j, a, b)
                exchange_integral = self.repulsion_dictionary.item(i, b, a, j)
                g_ij += self.density_matrix_total.item(a, b) * coulomb_integral - self.density_matrix.item(a, b) * exchange_integral
        return g_ij
