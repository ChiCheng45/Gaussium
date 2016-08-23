import itertools


class SinglesDoublesAmplitudes:

    def __init__(self, spin_molecular_integral, orbital_energies, occupied_orbitals, unoccupied_orbitals):
        self.spin_molecular_integral = spin_molecular_integral
        self.orbital_energies = orbital_energies
        self.occupied_orbitals = occupied_orbitals
        self.unoccupied_orbitals = unoccupied_orbitals
        self.total_orbitals = occupied_orbitals + unoccupied_orbitals
        self.singles, self.doubles, self.index_dict_1, self.index_dict_2 = self.indexes()

    def indexes(self):
        x = 0
        index_dict_1 = {}
        index_dict_2 = {}
        for i in range(self.occupied_orbitals):
            for a in range(self.occupied_orbitals, self.total_orbitals):
                index_dict_1[x] = (i, a)
                index_dict_2[(i, a)] = x
                x += 1
        singles = x
        for i, j in itertools.permutations(range(self.occupied_orbitals), 2):
            for a, b in itertools.permutations(range(self.occupied_orbitals, self.total_orbitals), 2):
                index_dict_1[x] = (i, j, a, b)
                index_dict_2[(i, j, a, b)] = x
                x += 1
        doubles = x - singles
        return singles, doubles, index_dict_1, index_dict_2

    def mp2_initial_guess(self):
        amplitudes = []
        for x in range(self.singles + self.doubles):

            if x <= self.singles - 1:
                amplitudes.append(0)
            else:
                i, j, a, b = self.index_dict_1[x]

                t_ijab = self.spin_molecular_integral.item(i, j, a, b) / (self.orbital_energies[i]
                + self.orbital_energies[j] - self.orbital_energies[a] - self.orbital_energies[b])

                amplitudes.append(t_ijab)

        return amplitudes

    def calculate_correlation(self, amplitudes):
        correlation = 0
        for i, j in itertools.permutations(range(self.occupied_orbitals), 2):
            for a, b in itertools.permutations(range(self.occupied_orbitals, self.total_orbitals), 2):

                index_ijab = self.index_dict_2[(i, j, a, b)]
                index_ia = self.index_dict_2[(i, a)]
                index_jb = self.index_dict_2[(j, b)]

                correlation += self.spin_molecular_integral.item(i, j, a, b) * (0.25 * amplitudes[index_ijab]
                + 0.5 * amplitudes[index_ia] * amplitudes[index_jb])

        return correlation
