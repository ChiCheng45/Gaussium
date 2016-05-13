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

    def create(self, orbital_coefficient, electrons):
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


class BlockedDensityMatrixUnrestricted:

    def __init__(self):
        self.matrix_size = 0
        self.half_matrix_size = 0
        self.electrons_alph = 0
        self.electrons_beta = 0
        self.orbital_coefficient = np.matrix([])

    def create(self, orbital_coefficient, electrons_alph, electrons_beta):
        self.electrons_alph = electrons_alph
        self.electrons_beta = electrons_beta
        self.orbital_coefficient = orbital_coefficient
        self.matrix_size = orbital_coefficient.shape[0]
        self.half_matrix_size = self.matrix_size // 2
        density_matrix = self.density_alph() + self.density_beta()
        return density_matrix

    def density_alph(self):
        density_matrix = np.zeros((self.matrix_size, self.matrix_size))
        for i in range(self.electrons_alph):
            density_matrix += self.orbital_coefficient[:, i] * self.orbital_coefficient[:, i].T
        return density_matrix

    def density_beta(self):
        density_matrix = np.zeros((self.matrix_size, self.matrix_size))
        for i in range(self.half_matrix_size, self.half_matrix_size + self.electrons_beta):
            density_matrix += self.orbital_coefficient[:, i] * self.orbital_coefficient[:, i].T
        return density_matrix
