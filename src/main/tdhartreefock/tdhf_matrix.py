from src.main.matrixelements import Matrix
from src.main.common import Indices


class TDHFMatrix(Matrix, Indices):

    def __init__(self, spin_molecular_integral, orbital_energies, occupied_orbitals, unoccupied_orbitals):
        Matrix.__init__(self, occupied_orbitals * unoccupied_orbitals)
        Indices.__init__(self, occupied_orbitals, unoccupied_orbitals)
        self.spin_molecular_integral = spin_molecular_integral
        self.orbital_energies = orbital_energies
        self.occupied_orbitals = occupied_orbitals
        self.unoccupied_orbitals = unoccupied_orbitals
        self.total_orbitals = occupied_orbitals + unoccupied_orbitals
        self.indices = self.pair_index()

    def pair_index(self):
        index = []
        for i, a in self.singles():
            index.append((i, a))
        return index

    def create_matrices(self):

        def calculate_a_elements(k, l):
            i, a = self.indices[k]
            j, b = self.indices[l]

            element = (i == j) * (a == b) * (self.orbital_energies[a] - self.orbital_energies[i])
            element += self.spin_molecular_integral[i, a, j, b] - self.spin_molecular_integral[i, j, a, b]
            return element

        def calculate_b_elements(k, l):
            i, a = self.indices[k]
            j, b = self.indices[l]

            element = self.spin_molecular_integral[i, a, j, b] - self.spin_molecular_integral[i, b, j, a]
            return element

        return self.create_matrix(calculate_a_elements), self.create_matrix(calculate_b_elements)


class TDHFMatrixSymmetryRestricted(Matrix, Indices):

    def __init__(self, spin_molecular_integral, orbital_energies, occupied_orbitals, unoccupied_orbitals):
        Matrix.__init__(self, occupied_orbitals * unoccupied_orbitals // 4)
        Indices.__init__(self, occupied_orbitals, unoccupied_orbitals)
        self.spin_molecular_integral = spin_molecular_integral
        self.orbital_energies = orbital_energies
        self.occupied_orbitals = occupied_orbitals
        self.unoccupied_orbitals = unoccupied_orbitals
        self.total_orbitals = occupied_orbitals + unoccupied_orbitals
        self.indices = self.pair_index()

    def pair_index(self):
        index = []
        for i, a in self.singles():
            if i % 2 == a % 2 == 0:
                index.append((i, a))
        return index

    def create_matrices_singlet(self):

        def calculate_a_elements(k, l):
            i, a = self.indices[k]
            j, b = self.indices[l]

            element = (i == j) * (a == b) * (self.orbital_energies[a] - self.orbital_energies[i])
            element += 2 * self.spin_molecular_integral[i, a, j, b] - self.spin_molecular_integral[i, j, a, b]
            return element

        def calculate_b_elements(k, l):
            i, a = self.indices[k]
            j, b = self.indices[l]

            element = 2 * self.spin_molecular_integral[i, a, j, b] - self.spin_molecular_integral[i, b, j, a]
            return element

        return self.create_matrix(calculate_a_elements), self.create_matrix(calculate_b_elements)

    def create_matrices_triplet(self):

        def calculate_a_elements(k, l):
            i, a = self.indices[k]
            j, b = self.indices[l]

            element = (i == j) * (a == b) * (self.orbital_energies[a] - self.orbital_energies[i])
            element -= self.spin_molecular_integral[i, j, a, b]
            return element

        def calculate_b_elements(k, l):
            i, a = self.indices[k]
            j, b = self.indices[l]

            element = - self.spin_molecular_integral[i, b, j, a]
            return element

        return self.create_matrix(calculate_a_elements), self.create_matrix(calculate_b_elements)
