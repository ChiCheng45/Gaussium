import numpy as np


class Matrix:

    def __init__(self, matrix_size):
        self.matrix_size = matrix_size

    def create_matrix(self, element):
        matrix = []
        for j in range(0, self.matrix_size):
            row = []
            for i in range(0, self.matrix_size):
                if j <= i:
                    matrix_ij = element.calculate(i, j)
                    row.append(matrix_ij)
                else:
                    row.append(0)
            matrix.append(row)
        matrix = np.matrix(matrix)
        matrix = matrix + matrix.T - np.diag(np.diag(matrix))
        return matrix
