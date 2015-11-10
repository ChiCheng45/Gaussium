import numpy as np
from scipy.spatial import distance
from src.main.common import VectorManipulation

class CoulombsLawArray:

    def __init__(self, nuclei_array):
        self.nuclei_array = nuclei_array

    def calculate_total_electric_potential_energy(self):
        energy_matrix = []
        for i in range(0, len(self.nuclei_array)):
            energy_matrix_row = []
            for j in range(0, len(self.nuclei_array)):
                if i == j:
                    energy_matrix_row.append(0)
                else:
                    energy = self.calc_electric_potential_energy(self.nuclei_array[i], self.nuclei_array[j])
                    energy_matrix_row.append(energy)
            energy_matrix.append(energy_matrix_row)
        return np.matrix(energy_matrix)

    def calc_electric_potential_energy(self, nuc1, nuc2):
        r_12 = VectorManipulation.squared_distance(nuc1.coordinates, nuc2.coordinates)
        electric_potential_energy = (nuc1.charge * nuc2.charge) / r_12
        return electric_potential_energy
