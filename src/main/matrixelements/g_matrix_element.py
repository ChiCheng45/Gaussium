class GMatrixElement:

    def __init__(self, density_matrix, repulsion_dictionary):
        self.density_matrix = density_matrix
        self.repulsion_dictionary = repulsion_dictionary

    def calculate(self, i, j):
        g_ij = 0
        for a in range(self.density_matrix.shape[0]):
            for b in range(self.density_matrix.shape[0]):
                coulomb_integral = self.repulsion_dictionary[self.sort((i, j, a, b))]
                exchange_integral = self.repulsion_dictionary[self.sort((i, b, a, j))]
                g_ij += self.density_matrix.item(b, a) * (coulomb_integral - (1/2) * exchange_integral)
        return g_ij

    @staticmethod
    def sort(i):
        a = i[0]
        b = i[1]
        c = i[2]
        d = i[3]
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
