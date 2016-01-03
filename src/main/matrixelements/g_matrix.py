from src.main.common import Symmetry
from src.main.matrixelements import Matrix
import numpy as np


class GMatrixRestricted(Matrix):

    def __init__(self, repulsion_dictionary):
        super().__init__()
        self.repulsion_dictionary = repulsion_dictionary
        self.density_matrix_total = np.matrix([])
        self.density_matrix = np.matrix([])

    def create_restricted(self, density_matrix_total):
        self.density_matrix_total = density_matrix_total
        self.matrix_size = density_matrix_total.shape[0]
        return self.create_matrix(self.calculate_restricted)

    def calculate_restricted(self, i, j):
        g_ij = 0
        for a in range(self.matrix_size):
            for b in range(self.matrix_size):
                coulomb_integral = self.repulsion_dictionary[Symmetry.sort(i, j, a, b)]
                exchange_integral = self.repulsion_dictionary[Symmetry.sort(i, b, a, j)]
                g_ij += self.density_matrix_total.item(b, a) * (coulomb_integral - (1/2) * exchange_integral)
        return g_ij


class GMatrixUnrestricted(Matrix):

    def __init__(self, repulsion_dictionary):
        super().__init__()
        self.repulsion_dictionary = repulsion_dictionary
        self.density_matrix_total = np.matrix([])
        self.density_matrix = np.matrix([])

    def create_unrestricted(self, density_matrix_1, density_matrix_2):
        self.density_matrix_total = density_matrix_1 + density_matrix_2
        self.density_matrix = density_matrix_1
        self.matrix_size = density_matrix_1.shape[0]
        return self.create_matrix(self.calculate_unrestricted)

    def calculate_unrestricted(self, i, j):
        g_ij = 0
        for a in range(self.matrix_size):
            for b in range(self.matrix_size):
                coulomb_integral = self.repulsion_dictionary[Symmetry.sort(i, j, a, b)]
                exchange_integral = self.repulsion_dictionary[Symmetry.sort(i, b, a, j)]
                g_ij += self.density_matrix_total.item(b, a) * coulomb_integral - self.density_matrix.item(b, a) * exchange_integral
        return g_ij
