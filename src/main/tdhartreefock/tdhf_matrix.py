from src.main.matrixelements import Matrix


class TDHFMatrix(Matrix):

    def __init__(self, spin_molecular_integral, orbital_energies, occupied_orbitals, unoccupied_orbitals):
        super().__init__(occupied_orbitals * unoccupied_orbitals)
        self.spin_molecular_integral = spin_molecular_integral
        self.orbital_energies = orbital_energies
        self.occupied_orbitals = occupied_orbitals
        self.unoccupied_orbitals = unoccupied_orbitals
        self.total_orbitals = occupied_orbitals + unoccupied_orbitals
        self.index_dict = self.pair_index()

    def pair_index(self):
        x = 0
        index_dict = {}
        for i in range(self.occupied_orbitals):
            for a in range(self.occupied_orbitals, self.total_orbitals):
                index_dict[x] = (i, a)
                x += 1
        return index_dict

    def create_matrices(self):

        def calculate_a_elements(k, l):
            i, a = self.index_dict[k]
            j, b = self.index_dict[l]

            element = (i == j) * (a == b) * (self.orbital_energies[a] - self.orbital_energies[i])
            element += self.spin_molecular_integral[i, a, j, b] - self.spin_molecular_integral[i, j, a, b]
            return element

        def calculate_b_elements(k, l):
            i, a = self.index_dict[k]
            j, b = self.index_dict[l]

            element = self.spin_molecular_integral[i, a, j, b] - self.spin_molecular_integral[i, b, j, a]
            return element

        return self.create_matrix(calculate_a_elements), self.create_matrix(calculate_b_elements)


class TDHFMatrixSymmetryRestricted(Matrix):

    def __init__(self, spin_molecular_integral, orbital_energies, occupied_orbitals, unoccupied_orbitals):
        super().__init__(occupied_orbitals * unoccupied_orbitals // 4)
        self.spin_molecular_integral = spin_molecular_integral
        self.orbital_energies = orbital_energies
        self.occupied_orbitals = occupied_orbitals
        self.unoccupied_orbitals = unoccupied_orbitals
        self.total_orbitals = occupied_orbitals + unoccupied_orbitals
        self.index_dict = self.pair_index()

    def pair_index(self):
        x = 0
        index_dict = {}
        for i in range(self.occupied_orbitals):
            for a in range(self.occupied_orbitals, self.total_orbitals):
                if i % 2 == a % 2 == 0:
                    index_dict[x] = (i, a)
                    x += 1
        return index_dict

    def create_matrices_singlet(self):

        def calculate_a_elements(k, l):
            i, a = self.index_dict[k]
            j, b = self.index_dict[l]

            element = (i == j) * (a == b) * (self.orbital_energies[a] - self.orbital_energies[i])
            element += 2 * self.spin_molecular_integral[i, a, j, b] - self.spin_molecular_integral[i, j, a, b]
            return element

        def calculate_b_elements(k, l):
            i, a = self.index_dict[k]
            j, b = self.index_dict[l]

            element = 2 * self.spin_molecular_integral[i, a, j, b] - self.spin_molecular_integral[i, b, j, a]
            return element

        return self.create_matrix(calculate_a_elements), self.create_matrix(calculate_b_elements)

    def create_matrices_triplet(self):

        def calculate_a_elements(k, l):
            i, a = self.index_dict[k]
            j, b = self.index_dict[l]

            element = (i == j) * (a == b) * (self.orbital_energies[a] - self.orbital_energies[i])
            element -= self.spin_molecular_integral[i, j, a, b]
            return element

        def calculate_b_elements(k, l):
            i, a = self.index_dict[k]
            j, b = self.index_dict[l]

            element = - self.spin_molecular_integral[i, b, j, a]
            return element

        return self.create_matrix(calculate_a_elements), self.create_matrix(calculate_b_elements)
