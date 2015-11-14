import numpy as np
from src.main.integrals.two_electron_repulsion_integral import ElectronRepulsionIntegral


class TwoElectronRepulsionElement:

    def __init__(self, basis_set_array):
        self.basis_set_array = basis_set_array

    def calculate(self, i, j, k, l):
        primitive_gaussian_array_i = self.basis_set_array[i].primitive_gaussian_array
        primitive_gaussian_array_j = self.basis_set_array[j].primitive_gaussian_array
        primitive_gaussian_array_k = self.basis_set_array[k].primitive_gaussian_array
        primitive_gaussian_array_l = self.basis_set_array[l].primitive_gaussian_array
        f_mn = 0
        for a in range(0, len(primitive_gaussian_array_i)):
            for b in range(0, len(primitive_gaussian_array_j)):
                for c in range(0, len(primitive_gaussian_array_k)):
                    for d in range(0, len(primitive_gaussian_array_l)):
                        a_1 = primitive_gaussian_array_i[a].exponent
                        a_2 = primitive_gaussian_array_j[b].exponent
                        a_3 = primitive_gaussian_array_k[c].exponent
                        a_4 = primitive_gaussian_array_l[d].exponent
                        c_1 = primitive_gaussian_array_i[a].contraction
                        c_2 = primitive_gaussian_array_j[b].contraction
                        c_3 = primitive_gaussian_array_k[c].contraction
                        c_4 = primitive_gaussian_array_l[d].contraction
                        l_1 = primitive_gaussian_array_i[a].integral_exponents
                        l_2 = primitive_gaussian_array_j[b].integral_exponents
                        l_3 = primitive_gaussian_array_k[c].integral_exponents
                        l_4 = primitive_gaussian_array_l[d].integral_exponents
                        n_1 = (((2 * a_1) / np.pi)**(3/4)) * (((((8 * a_1)**(l_1[0] + l_1[1] + l_1[2])) * np.math.factorial(l_1[0]) * np.math.factorial(l_1[1]) * np.math.factorial(l_1[2])) / (np.math.factorial(2 * l_1[0]) * np.math.factorial(2 * l_1[1]) * np.math.factorial(2 * l_1[2])))**(1/2))
                        n_2 = (((2 * a_2) / np.pi)**(3/4)) * (((((8 * a_2)**(l_2[0] + l_2[1] + l_2[2])) * np.math.factorial(l_2[0]) * np.math.factorial(l_2[1]) * np.math.factorial(l_2[2])) / (np.math.factorial(2 * l_2[0]) * np.math.factorial(2 * l_2[1]) * np.math.factorial(2 * l_2[2])))**(1/2))
                        n_3 = (((2 * a_3) / np.pi)**(3/4)) * (((((8 * a_3)**(l_3[0] + l_3[1] + l_3[2])) * np.math.factorial(l_3[0]) * np.math.factorial(l_3[1]) * np.math.factorial(l_3[2])) / (np.math.factorial(2 * l_3[0]) * np.math.factorial(2 * l_3[1]) * np.math.factorial(2 * l_3[2])))**(1/2))
                        n_4 = (((2 * a_4) / np.pi)**(3/4)) * (((((8 * a_4)**(l_4[0] + l_4[1] + l_4[2])) * np.math.factorial(l_4[0]) * np.math.factorial(l_4[1]) * np.math.factorial(l_4[2])) / (np.math.factorial(2 * l_4[0]) * np.math.factorial(2 * l_4[1]) * np.math.factorial(2 * l_4[2])))**(1/2))
                        f_mn += c_1 * c_2 * c_3 * c_4 * n_1 * n_2 * n_3 * n_4 * ElectronRepulsionIntegral.integral(primitive_gaussian_array_i[a], primitive_gaussian_array_j[b], primitive_gaussian_array_k[c], primitive_gaussian_array_l[d])
        return f_mn

    def store_integrals(self):
        electron_repulsion_dict = {}
        for a in range(0, len(self.basis_set_array)):
            for b in range(0, len(self.basis_set_array)):
                for c in range(0, len(self.basis_set_array)):
                    for d in range(0, len(self.basis_set_array)):
                        if not (a > b or c > d or a > c or (a == c and b > d)):
                            electron_repulsion_dict[(a, b, c, d)] = self.calculate(a, b, c, d)
        return electron_repulsion_dict
