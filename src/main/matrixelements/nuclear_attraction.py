import numpy as np
import scipy.special as sp
from src.main.common.vector_manipulation import VectorManipulation


class NuclearAttractionIntegral:

    def __init__(self, nuclei_array, basis_set_array):
        self.nuclei_array = nuclei_array
        self.basis_set_array = basis_set_array

    def calculate(self, i, j):
        v_ij = 0
        primitive_gaussian_array_i = self.basis_set_array[i].primitive_gaussian_array
        primitive_gaussian_array_j = self.basis_set_array[j].primitive_gaussian_array
        for a in range(0, len(primitive_gaussian_array_i)):
            for b in range(0, len(primitive_gaussian_array_j)):
                if primitive_gaussian_array_i[a].orbital_type == 'S' and primitive_gaussian_array_j[b].orbital_type == 'S':
                    a_1 = primitive_gaussian_array_i[a].exponent
                    a_2 = primitive_gaussian_array_j[b].exponent
                    d_1 = primitive_gaussian_array_i[a].contraction
                    d_2 = primitive_gaussian_array_j[b].contraction

                    r_1 = primitive_gaussian_array_i[a].coordinates
                    r_2 = primitive_gaussian_array_j[b].coordinates
                    r_12 = VectorManipulation.squared_distance(r_1, r_2)

                    n_1 = ((2 * a_1) / np.pi)**(3/4)
                    n_2 = ((2 * a_2) / np.pi)**(3/4)
                    s_ij = d_1 * d_2 * n_1 * n_2 * (np.pi / (a_1 + a_2))**(3/2) * np.exp(- a_1 * a_2 * r_12**2 / (a_1 + a_2))
                    for k in range(0, len(self.nuclei_array)):
                        c_k = self.nuclei_array[k].charge
                        r_3 = VectorManipulation.vector_gaussian(a_1, r_1, a_2, r_2)
                        r_k = self.nuclei_array[k].coordinates
                        r_3k = VectorManipulation.squared_distance(r_3, r_k)

                        if r_3k > 0:
                            f_ij = (1/2) * (np.pi / ((a_1 + a_2) * r_3k**2))**(1/2) * sp.erf(((a_1 + a_2) * r_3k**2)**(1/2))
                        else:
                            f_ij = 1
                        v_ij += - 2 * c_k * np.sqrt((a_1 + a_2) / np.pi) * s_ij * f_ij

                else:
                    v_ij += 0
        return v_ij

