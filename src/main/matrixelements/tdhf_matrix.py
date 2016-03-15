from src.main.matrixelements import Matrix
import numpy as np


class TDHFMatrix(Matrix):

    def __init__(self, molecular_integrals, orbital_energies, electrons):
        super().__init__()
        self.molecular_integrals = molecular_integrals
        self.orbital_energies = np.matrix(np.diag(orbital_energies))
        self.electrons = electrons
        self.matrix_size = orbital_energies.shape[0] * 2
        self.key = self.create_index_key()

    def create_index_key(self):
        key = {}
        i = -1
        for a in range(0, self.electrons):
            for b in range(0, self.electrons):
                i += 1
                key[i] = (a, b)
        return key

    def a_matrix(self):
        return self.create_matrix(self.calculate_a)

    def b_matrix(self):
        return self.create_matrix(self.calculate_b)

    def calculate_a(self, r, s):
        i, a = self.key[r]
        j, b = self.key[s]
        out1 = self.orbital_energies.item(i, j) - self.orbital_energies.item(a, b)
        out2 = self.molecular_integrals.item(i, a, j, b) - self.molecular_integrals.item(i, j, a, b)
        return out1 + out2

    def calculate_b(self, r, s):
        i, a = self.key[r]
        j, b = self.key[s]
        return self.molecular_integrals.item(i, a, j, b) - self.molecular_integrals.item(i, b, j, a)
