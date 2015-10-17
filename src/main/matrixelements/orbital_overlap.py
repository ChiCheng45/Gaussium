import numpy as np


class OverlapIntegral:

    def __init__(self, basis_set_array):
        self.basis_set_array = basis_set_array

    def calculate(self, i, j):
        s_ij = 0
        if i == j:
            return 1.00
        else:
            basis_coefficients_i = self.basis_set_array[i].get_array_of_coefficients()
            basis_coefficients_j = self.basis_set_array[j].get_array_of_coefficients()
            r_ab = np.linalg.norm(self.basis_set_array[i].get_coordinates() - self.basis_set_array[j].get_coordinates())
            for a in range(0, len(basis_coefficients_i)):
                for b in range(0, len(basis_coefficients_j)):
                    if self.basis_set_array[i].get_orbital_type() == 'S' and self.basis_set_array[j].get_orbital_type() == 'S':
                        a_1 = basis_coefficients_i[a][1]
                        a_2 = basis_coefficients_j[b][1]
                        c_1 = basis_coefficients_i[a][0]
                        c_2 = basis_coefficients_j[b][0]
                        n_1 = ((2 * a_1) / np.pi)**(3/4)
                        n_2 = ((2 * a_2) / np.pi)**(3/4)
                        s_ij += c_1 * c_2 * n_1 * n_2 * (np.pi / (a_1 + a_2))**(3/2) * np.exp(-a_1 * a_2 * r_ab**2 / (a_1 + a_2))
                    else:
                        s_ij += 0
            return s_ij
