import numpy as np
from src.main.common import Vector


class CoulombsLawMatrix:

    @classmethod
    def total_electric_potential_energy(cls, nuclei_array):
        matrix_length = len(nuclei_array)
        energy_matrix = np.matrix(np.zeros((matrix_length, matrix_length)))
        for i in range(matrix_length):
            for j in range(matrix_length):
                if i == j:
                    energy_matrix[i, j] = 0
                else:
                    energy_matrix[i, j] = cls.electric_potential_energy(nuclei_array[i], nuclei_array[j])
        return np.matrix(energy_matrix)

    @staticmethod
    def electric_potential_energy(nuc1, nuc2):
        r_12 = Vector.distance(nuc1.coordinates, nuc2.coordinates)
        electric_potential_energy = (nuc1.charge * nuc2.charge) / r_12
        return electric_potential_energy
