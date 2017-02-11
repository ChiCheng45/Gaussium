from src.main.common import read_basis_set_file


class NelderMead:

    def __init__(self, basis_file, energy_object):
        self.basis_file = basis_file
        self.energy_object = energy_object

    def optimize(self, nuclei_array):
        basis_set = read_basis_set_file(self.basis_file, nuclei_array)
        self.energy_object.calculate_energy(nuclei_array, basis_set)
