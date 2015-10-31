import numpy as np
import scipy.special as sp
from scipy.spatial import distance


class NuclearAttractionIntegral:

    def __init__(self, nuclei_array, basis_set_array):
        self.nuclei_array = nuclei_array
        self.basis_set_array = basis_set_array

    def calculate(self, i, j):
        v_ij = 0
        primitive_gaussian_array_i = self.basis_set_array[i].get_primitive_gaussian_array()
        primitive_gaussian_array_j = self.basis_set_array[j].get_primitive_gaussian_array()
        r_1 = self.basis_set_array[i].get_coordinates()
        r_2 = self.basis_set_array[j].get_coordinates()
        r_12 = distance.euclidean(r_1, r_2)
        for a in range(0, len(primitive_gaussian_array_i)):
            for b in range(0, len(primitive_gaussian_array_j)):
                if primitive_gaussian_array_i[a].get_orbital_type() == 'S' and primitive_gaussian_array_j[b].get_orbital_type() == 'S':
                    a_1 = primitive_gaussian_array_i[a].get_exponent()
                    a_2 = primitive_gaussian_array_j[b].get_exponent()
                    d_1 = primitive_gaussian_array_i[a].get_contraction()
                    d_2 = primitive_gaussian_array_j[b].get_contraction()
                    n_1 = ((2 * a_1) / np.pi)**(3/4)
                    n_2 = ((2 * a_2) / np.pi)**(3/4)
                    s_ij = d_1 * d_2 * n_1 * n_2 * (np.pi / (a_1 + a_2))**(3/2) * np.exp(- a_1 * a_2 * r_12**2 / (a_1 + a_2))
                    for k in range(0, len(self.nuclei_array)):
                        c_k = self.nuclei_array[k].get_charge()
                        r_3 = (a_1 * r_1 + a_2 * r_2) / (a_1 + a_2)
                        r_k = self.nuclei_array[k].get_coordinates()
                        r_3k = distance.euclidean(r_3, r_k)

                        if r_3k > 0:
                            f_ij = (1/2) * (np.pi / ((a_1 + a_2) * r_3k**2))**(1/2) * sp.erf(((a_1 + a_2) * r_3k**2)**(1/2))
                        else:
                            f_ij = 1
                        v_ij += - 2 * c_k * np.sqrt((a_1 + a_2) / np.pi) * s_ij * f_ij

                else:
                    v_ij += 0
        return v_ij

