class TwoElectronPartOfTheFockMatrixElements:

    def __init__(self, density_matrix, tei_matrix):
        self.density_matrix = density_matrix
        self.tei_matrix = tei_matrix

    def calculate(self, i, j):
        g_ij = 0
        p = self.density_matrix
        a = self.density_matrix.shape[0]
        for b in range(0, a):
            for c in range(0, a):

                    if i == 0 and j == 0:
                        k = 0
                    elif (i == 1 and j == 0) or (j == 1 and i == 0):
                        k = 1
                    else:
                        k = 2

                    if b == 0 and c == 0:
                        l = 0
                    elif (b == 1 and c == 0) or (c == 1 and b == 0):
                        l = 1
                    else:
                        l = 2

                    if i == 0 and c == 0:
                        m = 0
                    elif (i == 1 and c == 0) or (c == 1 and i == 0):
                        m = 1
                    else:
                        m = 2

                    if b == 0 and j == 0:
                        n = 0
                    elif (b == 1 and j == 0) or (j == 1 and b == 0):
                        n = 1
                    else:
                        n = 2

                    g_ij += p.item((b, c)) * ((self.tei_matrix.item((0, 1)).item((k, l))) - (1/2) * (self.tei_matrix.item((0, 1)).item(m, n)))
        return g_ij
