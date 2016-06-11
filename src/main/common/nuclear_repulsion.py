from src.main.common import Vector
import numpy as np


class CoulombsLawMatrix:

    @staticmethod
    def coulombs_law(nuc1, nuc2):
        r_12 = Vector.distance(nuc1.coordinates, nuc2.coordinates)
        ans = (nuc1.charge * nuc2.charge) / r_12
        return ans

    @classmethod
    def create(cls, nuclei_array):
        matrix_length = len(nuclei_array)
        matrix = np.matrix(np.zeros((matrix_length, matrix_length)))
        for i in range(matrix_length):
            for j in range(matrix_length):
                if i != j:
                    matrix[i, j] = cls.coulombs_law(nuclei_array[i], nuclei_array[j])
        return matrix
