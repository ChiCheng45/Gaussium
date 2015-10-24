import numpy as np
from src.main.common import SFunction
from src.main.common import BinomialCoefficientsFunction

class OverlapIntegral:

    def __init__(self, basis_set_array):
        self.basis_set_array = basis_set_array

    def calculate(self, i, j):
        s_ij = 0
        if i == j:
            return 1.00
        else:
            primitive_gaussian_array_i = self.basis_set_array[i].get_primitive_gaussian_array()
            primitive_gaussian_array_j = self.basis_set_array[j].get_primitive_gaussian_array()
            r_a = self.basis_set_array[i].get_coordinates()
            r_b = self.basis_set_array[j].get_coordinates()
            r_ab = np.linalg.norm(r_a - r_b)
            for a in range(0, len(primitive_gaussian_array_i)):
                for b in range(0, len(primitive_gaussian_array_j)):
                    a_1 = primitive_gaussian_array_i[a].get_exponent()
                    a_2 = primitive_gaussian_array_j[b].get_exponent()
                    c_1 = primitive_gaussian_array_i[a].get_contraction()
                    c_2 = primitive_gaussian_array_j[b].get_contraction()
                    l_1 = primitive_gaussian_array_i[a].get_integral_exponents()
                    l_2 = primitive_gaussian_array_i[a].get_integral_exponents()

                    n_1 = ((2 * a_1) / np.pi)**(3/4)
                    n_2 = ((2 * a_2) / np.pi)**(3/4)

                    r_p = (a_1 * r_a + a_2 * r_b) / (a_1 + a_2)
                    r_pa_xyz = r_a - r_p
                    r_pb_xyz = r_b - r_a

                    gamma = a_1 + a_2

                    s = SFunction(BinomialCoefficientsFunction)
                    s_x = s.calc(l_1[0], l_2[0], r_pa_xyz.item(0), r_pb_xyz.item(0), gamma)
                    s_y = s.calc(l_1[1], l_2[1], r_pa_xyz.item(1), r_pb_xyz.item(1), gamma)
                    s_z = s.calc(l_1[2], l_2[2], r_pa_xyz.item(2), r_pb_xyz.item(2), gamma)
                    s_ij += c_1 * c_2 * n_1 * n_2 * (np.pi / gamma)**(3/2) * np.exp(- a_1 * a_2 * r_ab**2 / gamma) * s_x * s_y * s_z
            return s_ij
