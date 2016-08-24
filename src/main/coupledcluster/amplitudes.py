import itertools


class SinglesDoublesAmplitudes:

    def __init__(self, spin_molecular_integral, orbital_energies, occupied_orbitals, unoccupied_orbitals):
        self.spin_molecular_integral = spin_molecular_integral
        self.orbital_energies = orbital_energies
        self.occupied_orbitals = range(occupied_orbitals)
        self.unoccupied_orbitals = range(occupied_orbitals, occupied_orbitals + unoccupied_orbitals)
        self.singles, self.doubles = self.indexes()

    def indexes(self):
        singles = []
        doubles = []
        for i in self.occupied_orbitals:
            for a in self.unoccupied_orbitals:
                singles.append((i, a))
        for i, j in itertools.permutations(self.occupied_orbitals, 2):
            for a, b in itertools.permutations(self.unoccupied_orbitals, 2):
                doubles.append((i, j, a, b))
        return singles, doubles

    def calculate_correlation(self, t):
        correlation = 0
        for i, j, a, b in self.doubles:
            correlation += self.spin_molecular_integral.item(i, j, a, b) * (0.25 * t[i, j, a, b]
            + 0.5 * t[i, a] * t[j, b])
        return correlation

    def mp2_initial_guess(self):
        t = {}
        for i, a in self.singles:
            t[i, a] = 0
        for i, j, a, b in self.doubles:
            t[i, j, a, b] = self.spin_molecular_integral.item(i, j, a, b) / (self.orbital_energies[i]
            + self.orbital_energies[j] - self.orbital_energies[a] - self.orbital_energies[b])
        return t

    def calculate_amplitudes(self, t_prev):
        t = {}
        inter = self.intermediates(t_prev)

        for i, a in self.singles:
            t[i, a] = 0
        for i, j, a, b in self.doubles:
            t[i, j, a, b] = self.spin_molecular_integral.item(i, j, a, b) / (self.orbital_energies[i]
            + self.orbital_energies[j] - self.orbital_energies[a] - self.orbital_energies[b])

        return t

    def denominator_arrays(self):
        pass

    def intermediates(self, t):
        intermediates = {}
        tau_1, tau_2 = self.tau(t)

        for a, e in itertools.permutations(self.unoccupied_orbitals, 2):
            f_ae = 0
            for m, f in self.singles:
                f_ae += t[m, f] * self.spin_molecular_integral.item(m, a, f, e)
                for n in self.occupied_orbitals:
                    f_ae -= 0.5 * tau_2[m, n, a, f] * self.spin_molecular_integral.item(m, n, e, f)
            intermediates[a, e] = f_ae

        for m, i in itertools.permutations(self.occupied_orbitals, 2):
            f_mi = 0
            for n, e in self.singles:
                f_mi += t[n, e] * self.spin_molecular_integral.item(m, n, i, e)
                for f in self.unoccupied_orbitals:
                    f_mi += 0.5 * tau_2[i, n, e, f] * self.spin_molecular_integral.item(m, n, e, f)
            intermediates[m, i] = f_mi

        for m, e in itertools.permutations(self.occupied_orbitals, self.unoccupied_orbitals):
            f_me = 0
            for n, f in self.singles:
                f_me += t[n, f] * self.spin_molecular_integral.item(m, n, e, f)
            intermediates[m, e] = f_me

        for m, n, i, j in itertools.permutations(self.occupied_orbitals, 4):
            w_mnij = self.spin_molecular_integral.item(m, n, i, j)
            for e in self.unoccupied_orbitals:
                w_mnij += t[j, e] * self.spin_molecular_integral.item(m, n, i, e) \
                - t[i, e] * self.spin_molecular_integral.item(m, n, j, e)
                for f in self.unoccupied_orbitals:
                    w_mnij += 0.25 * tau_1[i, j, e, f] * self.spin_molecular_integral.item(m, n, e, f)
            intermediates[m, n, i, j] = w_mnij

        for a, b, e, f in itertools.permutations(self.unoccupied_orbitals, 4):
            w_abef = self.spin_molecular_integral.item(a, b, e, f)
            for m in self.occupied_orbitals:
                w_abef -= t[m, b] * self.spin_molecular_integral.item(a, m, e, f) \
                - t[m, a] * self.spin_molecular_integral.item(b, m, e, f)
                for n in self.occupied_orbitals:
                    w_abef += 0.25 * tau_1[m, n, a, b] * self.spin_molecular_integral.item(m, n, e, f)
            intermediates[a, b, e, f] = w_abef

        for m, b in itertools.permutations(self.occupied_orbitals, self.unoccupied_orbitals):
            for e, j in itertools.permutations(self.unoccupied_orbitals, self.occupied_orbitals):
                w_mbej = self.spin_molecular_integral.item(m, b, e, j)
                for f in self.unoccupied_orbitals:
                    w_mbej += t[j, f] * self.spin_molecular_integral.item(m, b, e, f)
                for n in self.occupied_orbitals:
                    w_mbej -= t[n, b] * self.spin_molecular_integral.item(m, n, e, j)
                    for f in self.unoccupied_orbitals:
                        w_mbej -= (0.5 * t[j, n, f, b] + t[j, f] * t[n, b]) \
                        * self.spin_molecular_integral.item(m, n, e, f)
                intermediates[m, b, e, j] = w_mbej

        return intermediates

    def tau(self, t):
        tau_1 = {}
        tau_2 = {}
        for i, j, a, b in self.doubles:
            tau_1[i, j, a, b] = t[i, j, a, b] + t[i, a] * t[j, b] - t[i, b] * t[j, a]
            tau_2[i, j, a, b] = t[i, j, a, b] + 0.5 * (t[i, a] * t[j, b] - t[i, b] * t[j, a])
        return tau_1, tau_2
