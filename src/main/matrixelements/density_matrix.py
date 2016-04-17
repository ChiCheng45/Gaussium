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


class DensityMatrixRestrictedOpenShell(Matrix):

    def __init__(self):
        super().__init__()
        self.doubly_occupied_orbitals = 0
        self.singly_occupied_orbitals = 0
        self.unoccupied_orbitals = 0
        self.multiplicity = 0
        self.orbital_coefficient = np.matrix([])

    def create(self, orbital_coefficient, electrons, multiplicity):
        self.matrix_size = orbital_coefficient.shape[0]
        self.doubly_occupied_orbitals = (electrons - multiplicity + 1) // 2
        self.singly_occupied_orbitals = multiplicity - 1
        self.unoccupied_orbitals = self.matrix_size - self.doubly_occupied_orbitals - self.singly_occupied_orbitals
        self.orbital_coefficient = orbital_coefficient
        return self.create_matrix(self.calculate_doubly), self.create_matrix(self.calculate_singly), \
            self.create_matrix(self.calculate_unoccupied)

    def calculate_doubly(self, i, j):
        p_ij = 0
        for a in range(self.doubly_occupied_orbitals):
            p_ij += self.orbital_coefficient.item(i, a) * self.orbital_coefficient.item(j, a)
        return p_ij

    def calculate_singly(self, i, j):
        p_ij = 0
        start = self.doubly_occupied_orbitals
        end = self.doubly_occupied_orbitals + self.singly_occupied_orbitals
        for a in range(start, end):
            p_ij += self.orbital_coefficient.item(i, a) * self.orbital_coefficient.item(j, a)
        return p_ij

    def calculate_unoccupied(self, i, j):
        p_ij = 0
        start = self.singly_occupied_orbitals + self.doubly_occupied_orbitals
        end = self.matrix_size
        for a in range(start, end):
            p_ij += self.orbital_coefficient.item(i, a) * self.orbital_coefficient.item(j, a)
        return p_ij
