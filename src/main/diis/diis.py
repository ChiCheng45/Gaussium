import numpy as np


class DIIS:

    def __init__(self, overlap):
        self.overlap = overlap
        self.matrix_size = overlap.shape[0]
        self.fock_array = []
        self.error_array = []

    def fock_matrix(self, fock, density):
        error = self.overlap * density * fock - fock * density * self.overlap
        self.error_array.append(error)
        self.fock_array.append(fock)
        array_length = len(self.fock_array)

        if array_length > 1:
            diis_fock_matrix = np.matrix(np.zeros((self.matrix_size, self.matrix_size)))
            b_matrix = np.matrix(np.zeros((array_length + 1, array_length + 1)))
            vector_k = np.matrix(np.zeros(array_length + 1))
            vector_k.itemset(array_length, 1)

            for i in range(array_length):
                b_matrix.itemset((i, array_length), 1)
                b_matrix.itemset((array_length, i), 1)
                for j in range(array_length):
                    error_i = self.error_array[i]
                    error_j = self.error_array[j]
                    b_matrix.itemset((i, j), np.trace(error_i * np.transpose(error_j)))

            diis_coefficients = np.linalg.solve(b_matrix, np.transpose(vector_k))

            for l in range(array_length):
                diis_fock_matrix += diis_coefficients[l, 0] * self.fock_array[l]

            return diis_fock_matrix
        else:
            return fock
