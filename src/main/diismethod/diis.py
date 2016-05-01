import numpy as np


class DIIS:

    def __init__(self, overlap, linear_algebra):
        self.overlap = overlap
        self.linear_algebra = linear_algebra
        self.matrix_size = overlap.shape[0]
        self.fock_array = []
        self.error_array = []
        self.begin = False

    def fock_matrix(self, fock, density):
        error = self.overlap * density * fock - fock * density * self.overlap
        error = self.linear_algebra.orthonormalize(error)

        if np.any(error < 0.1 * np.ones((self.matrix_size, self.matrix_size))):
            self.begin = True
        if np.all(error < 1e-6 * np.ones((self.matrix_size, self.matrix_size))):
            self.begin = False

        if self.begin:
            self.error_array.append(error)
            self.fock_array.append(fock)

            if len(self.fock_array) <= 12:  # DIIS subspace reset

                if len(self.fock_array) > 1:
                    diis_fock_matrix = np.matrix(np.zeros((self.matrix_size, self.matrix_size)))
                    diis_coefficients = self.create_b_matrix()

                    for l in range(len(self.fock_array)):
                        diis_fock_matrix += diis_coefficients[l, 0] * self.fock_array[l]

                    return diis_fock_matrix
                else:
                    return fock

            else:
                self.fock_array = []
                self.error_array = []
                return fock

        else:
            return fock

    def create_b_matrix(self):
        array_length = len(self.fock_array)

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

        try:
            diis_coefficients = np.linalg.solve(b_matrix, np.transpose(vector_k))
            return diis_coefficients
        except np.linalg.linalg.LinAlgError:
            self.fock_array.pop(0)
            self.error_array.pop(0)
            print("np.linalg.linalg.LinAlgError: removing a DIIS vector")
            return self.create_b_matrix()
