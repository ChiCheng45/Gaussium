import numpy as np
from src.main.integrals import NuclearAttractionIntegral


class NuclearAttractionElement:

    def __init__(self, nuclei_array, basis_set_array):
        self.nuclei_array = nuclei_array
        self.basis_set_array = basis_set_array

    def calculate(self, i, j):
        v_ij = 0
        primitive_gaussian_array_i = self.basis_set_array[i].primitive_gaussian_array
        primitive_gaussian_array_j = self.basis_set_array[j].primitive_gaussian_array
        for a in range(len(primitive_gaussian_array_i)):
            for b in range(len(primitive_gaussian_array_j)):
                a_1 = primitive_gaussian_array_i[a].exponent
                a_2 = primitive_gaussian_array_j[b].exponent
                c_1 = primitive_gaussian_array_i[a].contraction
                c_2 = primitive_gaussian_array_j[b].contraction
                l_1 = primitive_gaussian_array_i[a].integral_exponents
                l_2 = primitive_gaussian_array_j[b].integral_exponents
                n_1 = (((2 * a_1) / np.pi)**(3/4)) * (((((8 * a_1)**(l_1[0] + l_1[1] + l_1[2])) * np.math.factorial(l_1[0]) * np.math.factorial(l_1[1]) * np.math.factorial(l_1[2])) / (np.math.factorial(2 * l_1[0]) * np.math.factorial(2 * l_1[1]) * np.math.factorial(2 * l_1[2])))**(1/2))
                n_2 = (((2 * a_2) / np.pi)**(3/4)) * (((((8 * a_2)**(l_2[0] + l_2[1] + l_2[2])) * np.math.factorial(l_2[0]) * np.math.factorial(l_2[1]) * np.math.factorial(l_2[2])) / (np.math.factorial(2 * l_2[0]) * np.math.factorial(2 * l_2[1]) * np.math.factorial(2 * l_2[2])))**(1/2))
                for k in range(len(self.nuclei_array)):
                    v_ij += - self.nuclei_array[k].charge * n_1 * n_2 * c_1 * c_2 * NuclearAttractionIntegral.primitive_nuclear_attraction(primitive_gaussian_array_i[a], primitive_gaussian_array_j[b], self.nuclei_array[k])
        return v_ij

