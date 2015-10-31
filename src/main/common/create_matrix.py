import numpy as np


class Matrix:

    def __init__(self, matrix_size):
        self.matrix_size = matrix_size

    def create_matrix(self, element):
        matrix = []
        for i in range(0, self.matrix_size):
            row = []
            for j in range(0, self.matrix_size):
                matrix_ij = element.calculate(i, j)
                row.append(matrix_ij)
            matrix.append(row)
        matrix = np.matrix(matrix)
        return matrix
