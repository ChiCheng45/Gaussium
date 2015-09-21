from src.main.Nuclei import Nuclei
import math


class Coulomb:
    def __init__(self, nuclei_array):
        self.nuclei_array = nuclei_array

    def calculate_total_electric_potential_energy(self):
        total_energy = 0
        num_nuclei = len(self.nuclei_array)
        for i in range(0, num_nuclei):
            for j in range(1, num_nuclei):
                if i < j:
                    energy = self.calculate_electric_potential_energy(self.nuclei_array[i], self.nuclei_array[j])
                    print(self.nuclei_array[i].get_name() + ' ' + self.nuclei_array[j].get_name() + ' ' + str(energy))
                    total_energy += energy
        return total_energy

    def calculate_electric_potential_energy(self, nuc_1: Nuclei, nuc_2: Nuclei):
        distance = math.sqrt((nuc_1.get_x() - nuc_2.get_x()) ** 2 + (nuc_1.get_y() - nuc_2.get_y()) ** 2 + (
            nuc_1.get_z() - nuc_2.get_z()) ** 2)
        electric_potential_energy = (nuc_1.get_charge() * nuc_2.get_charge()) / distance
        return electric_potential_energy
