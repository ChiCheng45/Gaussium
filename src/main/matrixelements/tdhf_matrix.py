from src.main.matrixelements import Matrix
import numpy as np


class TDHFMatrix(Matrix):

    def __init__(self, molecular_integrals, orbital_energies):
        super().__init__()
        self.molecular_integrals = molecular_integrals
        self.orbital_energies = np.matrix(np.diag(np.diag(orbital_energies)))
        self.matrix_size = molecular_integrals.shape[0]**2
        self.key = self.create_index_key()

    def create_index_key(self):
        integral_size = self.molecular_integrals.shape[0]
        key = {}
        for a in range(integral_size):
            for b in range(integral_size):
                for i in range(self.matrix_size):
                    key[i] = (a, b)
        return key

    def a_matrix(self):
        return self.create_matrix(self.calculate_a)

    def calculate_a(self, r, s):
        i, a = self.key[r]
        j, b = self.key[s]
        out1 = self.orbital_energies.item(i, j) - self.orbital_energies.item(a, b)
        out2 = self.molecular_integrals.item(i, a, j, b) - self.molecular_integrals.item(i, j, a, b)
        return out1 + out2

    def b_matrix(self):
        return self.create_matrix(self.calculate_b)

    def calculate_b(self, r, s):
        i, a = self.key[r]
        j, b = self.key[s]
        return self.molecular_integrals.item(i, a, j, b) - self.molecular_integrals.item(i, b, j, a)
