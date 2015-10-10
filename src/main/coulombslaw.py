import numpy as np


class CoulombTotal:

    def __init__(self, coulomb, nuclei_array):
        self.nuclei_array = nuclei_array
        self.coulomb = coulomb

    def calculate_total_electric_potential_energy(self):
        energy_matrix = []
        for i in range(0, len(self.nuclei_array)):
            energy_matrix_row = []
            for j in range(0, len(self.nuclei_array)):
                if i == j:
                    energy_matrix_row.append(0)
                else:
                    energy = self.coulomb(self.nuclei_array[i], self.nuclei_array[j]).calc_electric_potential_energy()
                    energy_matrix_row.append(energy)
            energy_matrix.append(energy_matrix_row)
        return np.matrix(energy_matrix)


class Coulomb:

    def __init__(self, nuc1, nuc2):
        self.nuc1 = nuc1
        self.nuc2 = nuc2

    def calc_electric_potential_energy(self):
        n1 = self.nuc1
        n2 = self.nuc2
        distance = np.sqrt((n1.get_x() - n2.get_x()) ** 2 + (n1.get_y() - n2.get_y()) ** 2 + (
            n1.get_z() - n2.get_z()) ** 2)
        electric_potential_energy = (n1.get_charge() * n2.get_charge()) / distance
        return electric_potential_energy
