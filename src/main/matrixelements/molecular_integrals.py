import numpy as np


class MolecularIntegrals:

    @staticmethod
    def calculate(repulsion, coefficients):
        matrix_size = coefficients.shape[0]
        for r in range(matrix_size):
            for s in range(matrix_size):
                repulsion[r, s, :, :] = np.transpose(coefficients) * repulsion[r, s, :, :] * coefficients
        for t in range(matrix_size):
            for u in range(matrix_size):
                repulsion[:, :, t, u] = np.transpose(coefficients) * repulsion[:, :, t, u] * coefficients
        return repulsion
