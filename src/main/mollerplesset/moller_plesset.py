class MollerPlesset:

    def __init__(self):
        pass

    @classmethod
    def second_order(self, electrons, orbital_energies, orbital_coefficients, repulsion):
        occupied_orbitals = electrons // 2
        ans = 0
        for i in range(occupied_orbitals):
            for j in range(occupied_orbitals):
                for a in range(occupied_orbitals, len(orbital_energies)):
                    for b in range(occupied_orbitals, len(orbital_energies)):
                        ans += (repulsion[self.sort(i, j, a, b)] - repulsion[self.sort(i, b, a, j)])**2 / (2*(orbital_energies[0] - orbital_energies[1]))
        return ans

    @staticmethod
    def sort(i, j, k, l):
        a = i
        b = j
        c = k
        d = l
        if a > b:
            a, b = b, a
        if c > d:
            c, d = d, c
        if a > c:
            a, c = c, a
            b, d = d, b
        if a == c:
            if b > d:
                a, c = c, a
                b, d = d, b
        return a, b, c, d