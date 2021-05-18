import itertools
from gaussium.common import Indices


class Amplitudes(Indices):

    def __init__(self, spin_molecular_integral, orbital_energies, occupied_orbitals, unoccupied_orbitals):
        super().__init__(occupied_orbitals, unoccupied_orbitals)
        self.integrals = spin_molecular_integral
        self.orbital_energies = orbital_energies
        self.denominator = {}

    def denominator_singles(self):
        for i, a in self.singles():
            self.denominator[i, a] = self.orbital_energies[i] - self.orbital_energies[a]

    def denominator_doubles(self):
        for i, j, a, b in self.doubles():
            self.denominator[i, j, a, b] = self.orbital_energies[i] + self.orbital_energies[j] \
            - self.orbital_energies[a] - self.orbital_energies[b]

    def denominator_triples(self):
        for i, j, k, a, b, c in self.triples():
            self.denominator[i, j, k, a, b, c] = self.orbital_energies[i] + self.orbital_energies[j] \
            + self.orbital_energies[k] - self.orbital_energies[a] - self.orbital_energies[b] - self.orbital_energies[c]


class SinglesDoubles(Amplitudes):

    def __init__(self, spin_molecular_integral, orbital_energies, occupied_orbitals, unoccupied_orbitals):
        super().__init__(spin_molecular_integral, orbital_energies, occupied_orbitals, unoccupied_orbitals)
        self.denominator_singles()
        self.denominator_doubles()

    def mp2_initial_guess(self):
        t = {}
        for i, a in self.singles():
            t[i, a] = 0
        for i, j, a, b in self.doubles():
            t[i, j, a, b] = self.integrals.item(i, j, a, b) / self.denominator[i, j, a, b]
        return t

    def calculate_amplitudes(self, t_prev):
        t = {}
        tau_1, tau_2 = self.tau(t_prev)
        inter = self.intermediates(t_prev, tau_1, tau_2)

        for i, a in self.singles():
            t[i, a] = self.singles_amplitudes(i, a, t_prev, inter)
        for i, j, a, b in self.doubles():
            t[i, j, a, b] = self.doubles_amplitudes(i, j, a, b, t_prev, tau_1, inter)

        return t

    def singles_amplitudes(self, i, a, t, inter):
        t_ia = 0

        for e in self.unoccupied_indices:
            t_ia += t[i, e] * inter[a, e]

        for m in self.occupied_indices:
            t_ia -= t[m, a] * inter[m, i]
            for e in self.unoccupied_indices:
                if i != m and a != e:
                    t_ia += t[i, m, a, e] * inter[m, e]

        for n, f in self.singles():
            t_ia -= t[n, f] * self.integrals.item(n, a, i, f)

        for m, e in self.singles():
            for f in self.unoccupied_indices:
                if i != m and e != f:
                    t_ia -= 0.5 * t[i, m, e, f] * self.integrals.item(m, a, e, f)
            for n in self.occupied_indices:
                if m != n and a != e:
                    t_ia -= 0.5 * t[m, n, a, e] * self.integrals.item(n, m, e, i)

        return t_ia / self.denominator[i, a]

    def doubles_amplitudes(self, i, j, a, b, t, tau, inter):
        t_ijab = self.integrals.item(i, j, a, b)

        for e in self.unoccupied_indices:
            out = 0
            for m in self.occupied_indices:
                out += t[m, b] * inter[m, e]
            if a != e:
                t_ijab += t[i, j, a, e] * (inter[b, e] - 0.5 * out)
            out = 0
            for m in self.occupied_indices:
                out += t[m, a] * inter[m, e]
            if b != e:
                t_ijab -= t[i, j, b, e] * (inter[a, e] - 0.5 * out)

        for m in self.occupied_indices:
            out = 0
            for e in self.unoccupied_indices:
                out += t[j, e] * inter[m, e]
            if i != m:
                t_ijab -= t[i, m, a, b] * (inter[m, j] + 0.5 * out)
            out = 0
            for e in self.unoccupied_indices:
                out += t[i, e] * inter[m, e]
            if j != m:
                t_ijab += t[j, m, a, b] * (inter[m, i] + 0.5 * out)

        for m, n in itertools.permutations(self.occupied_indices, 2):
            t_ijab += 0.5 * tau[m, n, a, b] * inter[m, n, i, j]

        for e, f in itertools.permutations(self.unoccupied_indices, 2):
            t_ijab += 0.5 * tau[i, j, e, f] * inter[a, b, e, f]

        for m, e in self.singles():
            if i != m and a != e:
                t_ijab += t[i, m, a, e] * inter[m, b, e, j]
            t_ijab -= t[i, e] * t[m, a] * self.integrals.item(m, b, e, j)
            if i != m and b != e:
                t_ijab -= t[i, m, b, e] * inter[m, a, e, j]
            t_ijab += t[i, e] * t[m, b] * self.integrals.item(m, a, e, j)
            if j != m and a != e:
                t_ijab -= t[j, m, a, e] * inter[m, b, e, i]
            t_ijab += t[j, e] * t[m, a] * self.integrals.item(m, b, e, i)
            if j != m and b != e:
                t_ijab += t[j, m, b, e] * inter[m, a, e, i]
            t_ijab -= t[j, e] * t[m, b] * self.integrals.item(m, a, e, i)

        for e in self.unoccupied_indices:
            t_ijab += t[i, e] * self.integrals.item(a, b, e, j) - t[j, e] * self.integrals.item(a, b, e, i)

        for m in self.occupied_indices:
            t_ijab -= t[m, a] * self.integrals.item(m, b, i, j) - t[m, b] * self.integrals.item(m, a, i, j)

        return t_ijab / self.denominator[i, j, a, b]

    def intermediates(self, t, tau_1, tau_2):
        intermediates = {}

        for a, e in itertools.product(self.unoccupied_indices, repeat=2):
            f_ae = 0
            for m, f in self.singles():
                f_ae += t[m, f] * self.integrals.item(m, a, f, e)
                for n in self.occupied_indices:
                    if m != n and a != f:
                        f_ae -= 0.5 * tau_2[m, n, a, f] * self.integrals.item(m, n, e, f)
            intermediates[a, e] = f_ae

        for m, i in itertools.product(self.occupied_indices, repeat=2):
            f_mi = 0
            for n, e in self.singles():
                f_mi += t[n, e] * self.integrals.item(m, n, i, e)
                for f in self.unoccupied_indices:
                    if i != n and e != f:
                        f_mi += 0.5 * tau_2[i, n, e, f] * self.integrals.item(m, n, e, f)
            intermediates[m, i] = f_mi

        for m, e in itertools.product(self.occupied_indices, self.unoccupied_indices):
            f_me = 0
            for n, f in self.singles():
                f_me += t[n, f] * self.integrals.item(m, n, e, f)
            intermediates[m, e] = f_me

        for m, n, i, j in itertools.product(self.occupied_indices, repeat=4):
            w_mnij = self.integrals.item(m, n, i, j)
            for e in self.unoccupied_indices:
                w_mnij += t[j, e] * self.integrals.item(m, n, i, e) - t[i, e] * self.integrals.item(m, n, j, e)
                for f in self.unoccupied_indices:
                    if i != j and e != f:
                        w_mnij += 0.25 * tau_1[i, j, e, f] * self.integrals.item(m, n, e, f)
            intermediates[m, n, i, j] = w_mnij

        for a, b, e, f in itertools.product(self.unoccupied_indices, repeat=4):
            w_abef = self.integrals.item(a, b, e, f)
            for m in self.occupied_indices:
                w_abef -= t[m, b] * self.integrals.item(a, m, e, f) - t[m, a] * self.integrals.item(b, m, e, f)
                for n in self.occupied_indices:
                    if m != n and a != b:
                        w_abef += 0.25 * tau_1[m, n, a, b] * self.integrals.item(m, n, e, f)
            intermediates[a, b, e, f] = w_abef

        for m, b, e, j in itertools.product(self.occupied_indices, self.unoccupied_indices, self.unoccupied_indices,
        self.occupied_indices):
            w_mbej = self.integrals.item(m, b, e, j)
            for f in self.unoccupied_indices:
                w_mbej += t[j, f] * self.integrals.item(m, b, e, f)
            for n in self.occupied_indices:
                w_mbej -= t[n, b] * self.integrals.item(m, n, e, j)
                for f in self.unoccupied_indices:
                    if j != n and f != b:
                        w_mbej -= 0.5 * t[j, n, f, b] * self.integrals.item(m, n, e, f)
                    w_mbej -= t[j, f] * t[n, b] * self.integrals.item(m, n, e, f)
            intermediates[m, b, e, j] = w_mbej

        return intermediates

    def tau(self, t):
        tau_1 = {}
        tau_2 = {}
        for i, j, a, b in self.doubles():
            tau_1[i, j, a, b] = t[i, j, a, b] + t[i, a] * t[j, b] - t[i, b] * t[j, a]
            tau_2[i, j, a, b] = t[i, j, a, b] + 0.5 * (t[i, a] * t[j, b] - t[i, b] * t[j, a])
        return tau_1, tau_2


class PeturbativeTriples(Amplitudes):

    def __init__(self, spin_molecular_integral, orbital_energies, occupied_orbitals, unoccupied_orbitals):
        super().__init__(spin_molecular_integral, orbital_energies, occupied_orbitals, unoccupied_orbitals)
        self.denominator_triples()

    def calculate_triples_amplitudes(self, t):
        t_connected = {}
        t_disconnected = {}
        for i, j, k, a, b, c in self.triples():
            t_connected[i, j, k, a, b, c] = self.connected_triples_amplitude(i, j, k, a, b, c, t)
            t_disconnected[i, j, k, a, b, c] = self.disconnected_triples_amplitude(i, j, k, a, b, c, t)
        return t_connected, t_disconnected

    def connected_triples_amplitude(self, i, j, k, a, b, c, t):
        numerator = self.connected_triples_inter(i, j, k, a, b, c, t) - self.connected_triples_inter(j, i, k, a, b, c, t) - self.connected_triples_inter(k, j, i, a, b, c, t) \
        - self.connected_triples_inter(i, j, k, b, a, c, t) + self.connected_triples_inter(j, i, k, b, a, c, t) + self.connected_triples_inter(k, j, i, b, a, c, t) \
        - self.connected_triples_inter(i, j, k, c, b, a, t) + self.connected_triples_inter(j, i, k, c, b, a, t) + self.connected_triples_inter(k, j, i, c, b, a, t)
        return numerator / self.denominator[i, j, k, a, b, c]

    def connected_triples_inter(self, i, j, k, a, b, c, t):
        ans = 0
        for e in self.unoccupied_indices:
            if a != e:
                ans += t[j, k, a, e] * self.integrals[e, i, b, c]
        for m in self.occupied_indices:
            if i != m:
                ans -= t[i, m, b, c] * self.integrals[m, a, j, k]
        return ans

    def disconnected_triples_amplitude(self, i, j, k, a, b, c, t):
        numerator = t[i, a] * self.integrals[j, k, b, c] - t[j, a] * self.integrals[i, k, b, c] - t[k, a] * self.integrals[j, i, b, c] \
        - t[i, b] * self.integrals[j, k, a, c] + t[j, b] * self.integrals[i, k, a, c] + t[k, b] * self.integrals[j, i, a, c] \
        - t[i, c] * self.integrals[j, k, b, a] + t[j, c] * self.integrals[i, k, b, a] + t[k, c] * self.integrals[j, i, b, a]
        return numerator / self.denominator[i, j, k, a, b, c]
