from src.main.matrixelements import Matrix
import numpy as np


class DensityMatrixRestricted(Matrix):

    def __init__(self, electrons):
        super().__init__()
        self.electrons = electrons
        self.orbital_coefficient = np.matrix([])

    def create(self, orbital_coefficient):
        self.orbital_coefficient = orbital_coefficient
        self.matrix_size = orbital_coefficient.shape[0]
        return self.create_matrix(self.calculate_restricted)

    def calculate_restricted(self, i, j):
        p_ij = 0
        c = self.orbital_coefficient
        for a in range(self.electrons // 2):
            p_ij += 2 * c.item(i, a) * c.item(j, a)
        return p_ij


class DensityMatrixUnrestricted(Matrix):

    def __init__(self):
        super().__init__()
        self.electrons = 0
        self.orbital_coefficient = np.matrix([])

    def create(self, electrons, orbital_coefficient):
        self.electrons = electrons
        self.orbital_coefficient = orbital_coefficient
        self.matrix_size = orbital_coefficient.shape[0]
        return self.create_matrix(self.calculate_unrestricted)

    def calculate_unrestricted(self, i, j):
        p_ij = 0
        c = self.orbital_coefficient
        for a in range(self.electrons):
            p_ij += c.item(i, a) * c.item(j, a)
        return p_ij


