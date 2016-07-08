import numpy as np


class Matrix:

    def __init__(self, matrix_size):
        self.matrix_size = matrix_size

    def create_matrix(self, method):
        matrix = np.matrix(np.zeros((self.matrix_size, self.matrix_size)))
        for i in range(self.matrix_size):
            for j in range(self.matrix_size):
                if i <= j:
                    matrix.itemset((i, j), method(i, j))
        return matrix + np.transpose(matrix) - np.diag(np.diag(matrix))
