class MolecularIntegrals:

    def __init__(self, repulsion, orbital_coefficients):
        self.repulsion = repulsion
        self.orbital_coefficients = orbital_coefficients

    def calc(self, i, j, k, l):
        ans = 0
        c = self.orbital_coefficients
        for r in range(c.shape[0]):
            for s in range(c.shape[0]):
                for t in range(c.shape[0]):
                    for u in range(c.shape[0]):
                        coefficients = c.item(r, i) * c.item(s, j) * c.item(t, k) * c.item(u, l)
                        ans += coefficients * self.repulsion[r, s, t, u]
        return ans
