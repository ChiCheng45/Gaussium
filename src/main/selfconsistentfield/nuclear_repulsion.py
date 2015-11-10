import numpy as np
from scipy.spatial import distance
from src.main.common import Vector


class CoulombsLawArray:

    @classmethod
    def calculate_total_electric_potential_energy(cls, nuclei_array):
        energy_matrix = []
        for i in range(0, len(nuclei_array)):
            energy_matrix_row = []
            for j in range(0, len(nuclei_array)):
                if i == j:
                    energy_matrix_row.append(0)
                else:
                    energy = cls.calc_electric_potential_energy(nuclei_array[i], nuclei_array[j])
                    energy_matrix_row.append(energy)
            energy_matrix.append(energy_matrix_row)
        return np.matrix(energy_matrix)

    @staticmethod
    def calc_electric_potential_energy(nuc1, nuc2):
        r_12 = Vector.distance(nuc1.coordinates, nuc2.coordinates)
        electric_potential_energy = (nuc1.charge * nuc2.charge) / r_12
        return electric_potential_energy
