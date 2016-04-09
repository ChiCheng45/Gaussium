from src.main.matrixelements import Matrix
import numpy as np


class GMatrixRestricted(Matrix):

    def __init__(self, repulsion_matrix):
        super().__init__()
        self.repulsion_matrix = repulsion_matrix
        self.density_matrix_total = np.matrix([])

    def create(self, density_matrix):
        self.matrix_size = density_matrix.shape[0]
        self.density_matrix_total = density_matrix
        return self.create_matrix(self.calculate_restricted)

    def calculate_restricted(self, i, j):
        g_ij = 0
        for a in range(self.matrix_size):
            for b in range(self.matrix_size):
                coulomb_integral = self.repulsion_matrix.item(i, j, a, b)
                exchange_integral = self.repulsion_matrix.item(i, b, a, j)
                g_ij += self.density_matrix_total.item(a, b) * (coulomb_integral - 1/2 * exchange_integral)
        return g_ij


class GMatrixUnrestricted(Matrix):

    def __init__(self, repulsion_matrix):
        super().__init__()
        self.repulsion_matrix = repulsion_matrix
        self.density_matrix_alph = np.matrix([])
        self.density_matrix_beta = np.matrix([])
        self.density_matrix_total = np.matrix([])

    def create(self, density_matrix_alph, density_matrix_beta):
        self.matrix_size = density_matrix_alph.shape[0]
        self.density_matrix_alph = density_matrix_alph
        self.density_matrix_beta = density_matrix_beta
        self.density_matrix_total = density_matrix_alph + density_matrix_beta
        return self.create_matrix(self.calc_unrestricted_alph), self.create_matrix(self.calc_unrestricted_beta)

    def calc_unrestricted_alph(self, i, j):
        g_ij = 0
        for a in range(self.matrix_size):
            for b in range(self.matrix_size):
                coulomb = self.density_matrix_total.item(a, b) * self.repulsion_matrix.item(i, j, a, b)
                exchange = self.density_matrix_alph.item(a, b) * self.repulsion_matrix.item(i, b, a, j)
                g_ij += coulomb - exchange
        return g_ij

    def calc_unrestricted_beta(self, i, j):
        g_ij = 0
        for a in range(self.matrix_size):
            for b in range(self.matrix_size):
                coulomb = self.density_matrix_total.item(a, b) * self.repulsion_matrix.item(i, j, a, b)
                exchange = self.density_matrix_beta.item(a, b) * self.repulsion_matrix.item(i, b, a, j)
                g_ij += coulomb - exchange
        return g_ij


class GMatrixConstrainedUnrestricted(Matrix):

    def __init__(self, repulsion_matrix, electrons_alph, electrons_beta):
        super().__init__()
        self.repulsion_matrix = repulsion_matrix
        self.density_matrix_alph = np.matrix([])
        self.density_matrix_beta = np.matrix([])
        self.density_matrix_total = np.matrix([])
        self.electrons_alph = electrons_alph
        self.electrons_beta = electrons_beta

    def create(self, density_matrix_alph, density_matrix_beta):
        self.matrix_size = density_matrix_alph.shape[0]
        self.density_matrix_alph = density_matrix_alph
        self.density_matrix_beta = density_matrix_beta
        self.density_matrix_total = density_matrix_alph + density_matrix_beta
        return self.create_matrix(self.calc_constrained_alph), self.create_matrix(self.calc_constrained_beta)

    def calc_constrained_alph(self, i, j):
        # (i is virtual and j is open or closed) or (i is virtual and j is open or closed)
        if i > self.electrons_alph and j < (self.electrons_alph or self.electrons_beta) \
        or j > self.electrons_alph and i < (self.electrons_alph or self.electrons_beta):

            g_ij = 0
            for a in range(self.matrix_size):
                for b in range(self.matrix_size):
                    coulomb = self.repulsion_matrix.item(i, j, a, b)
                    exchange = self.repulsion_matrix.item(i, b, a, j)
                    g_ij += self.density_matrix_total.item(a, b) * (coulomb - 1/2 * exchange)
            return g_ij

        else:

            g_ij = 0
            for a in range(self.matrix_size):
                for b in range(self.matrix_size):
                    coulomb = self.density_matrix_total.item(a, b) * self.repulsion_matrix.item(i, j, a, b)
                    exchange = self.density_matrix_alph.item(a, b) * self.repulsion_matrix.item(i, b, a, j)
                    g_ij += coulomb - exchange
            return g_ij

    def calc_constrained_beta(self, i, j):
        # (i is virtual and j is open or closed) or (i is virtual and j is open or closed)
        if i > self.electrons_beta and j < (self.electrons_alph or self.electrons_beta) \
        or j > self.electrons_beta and i < (self.electrons_alph or self.electrons_beta):

            g_ij = 0
            for a in range(self.matrix_size):
                for b in range(self.matrix_size):
                    coulomb = self.repulsion_matrix.item(i, j, a, b)
                    exchange = self.repulsion_matrix.item(i, b, a, j)
                    g_ij += self.density_matrix_total.item(a, b) * (coulomb - 1/2 * exchange)
            return g_ij

        else:

            g_ij = 0
            for a in range(self.matrix_size):
                for b in range(self.matrix_size):
                    coulomb = self.density_matrix_total.item(a, b) * self.repulsion_matrix.item(i, j, a, b)
                    exchange = self.density_matrix_beta.item(a, b) * self.repulsion_matrix.item(i, b, a, j)
                    g_ij += coulomb - exchange
            return g_ij
