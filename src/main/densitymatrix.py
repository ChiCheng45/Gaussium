class DensityMatrix:

    def __init__(self, orbital_coefficient_matrix):
        self.orbital_coefficient_matrix = orbital_coefficient_matrix

    def calculate(self, i, j):
        p_ij = 0
        c = self.orbital_coefficient_matrix
        electrons = 2
        for a in range(0, int((1/2) * electrons)):
            p_ij += 2 * c.item((i, a)) * c.item((j, a))
        return p_ij
