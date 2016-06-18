from src.main.common import coordinate_distance
import numpy as np


def coulomb_matrix(nuclei_array):
    matrix_length = len(nuclei_array)
    matrix = np.matrix(np.zeros((matrix_length, matrix_length)))
    for i in range(matrix_length):
        for j in range(matrix_length):
            if i != j:
                matrix[i, j] = coulombs_law(nuclei_array[i], nuclei_array[j])
    return matrix


def coulombs_law(nuc1, nuc2):
    r_12 = coordinate_distance(nuc1.coordinates, nuc2.coordinates)
    ans = (nuc1.charge * nuc2.charge) / r_12
    return ans
