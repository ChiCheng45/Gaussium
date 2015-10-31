import numpy as np


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
        distance = np.linalg.norm(nuc1.get_coordinates() - nuc2.get_coordinates())
        electric_potential_energy = (nuc1.get_charge() * nuc2.get_charge()) / distance
        return electric_potential_energy
