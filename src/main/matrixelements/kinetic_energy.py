import numpy as np


class KineticEnergyIntegral:

    def __init__(self, basis_set_array):
        self.basis_set_array = basis_set_array

    def calculate(self, i, j):
        t_ij = 0
        primitive_gaussian_array_i = self.basis_set_array[i].get_primitive_gaussian_array()
        primitive_gaussian_array_j = self.basis_set_array[j].get_primitive_gaussian_array()
        r_ab = np.linalg.norm(self.basis_set_array[i].get_coordinates() - self.basis_set_array[j].get_coordinates())
        for a in range(0, len(primitive_gaussian_array_i)):
            for b in range(0, len(primitive_gaussian_array_j)):
                if primitive_gaussian_array_i[a].get_orbital_type() == 'S' and primitive_gaussian_array_j[b].get_orbital_type() == 'S':
                    a_1 = primitive_gaussian_array_i[a].get_exponent()
                    a_2 = primitive_gaussian_array_j[b].get_exponent()
                    c_1 = primitive_gaussian_array_i[a].get_contraction()
                    c_2 = primitive_gaussian_array_j[b].get_contraction()
                    n_1 = ((2 * a_1) / np.pi)**(3/4)
                    n_2 = ((2 * a_2) / np.pi)**(3/4)
                    s_ij = c_1 * c_2 * n_1 * n_2 * (np.pi / (a_1 + a_2))**(3/2) * np.exp(- a_1 * a_2 * r_ab**2 / (a_1 + a_2))
                    k_ij = 3 * a_1 * a_2 / (a_1 + a_2) - 2 * a_1**2 * a_2**2 * r_ab**2 / (a_1 + a_2)**2
                    t_ij += s_ij * k_ij
                else:
                    t_ij += 0
        return t_ij
