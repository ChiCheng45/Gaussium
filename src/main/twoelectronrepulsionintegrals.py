import numpy as np
import scipy.special as sp

# some really bad code below, sorry

class TwoElectronRepulsion:

    def __init__(self, basis_set_array):
        self.basis_set_array = basis_set_array

    def calculate(self, i, j, k, l):
        basis_coefficients_i = self.basis_set_array[i].get_array_of_coefficients()
        basis_coefficients_j = self.basis_set_array[j].get_array_of_coefficients()
        basis_coefficients_k = self.basis_set_array[k].get_array_of_coefficients()
        basis_coefficients_l = self.basis_set_array[l].get_array_of_coefficients()
        f_mn = 0
        for a in range(0, len(basis_coefficients_i)):
            for b in range(0, len(basis_coefficients_j)):
                for c in range(0, len(basis_coefficients_k)):
                    for d in range(0, len(basis_coefficients_l)):

                        a_1 = float(basis_coefficients_i[a][1])
                        a_2 = float(basis_coefficients_j[b][1])
                        a_3 = float(basis_coefficients_k[c][1])
                        a_4 = float(basis_coefficients_l[d][1])

                        c_1 = float(basis_coefficients_i[a][0])
                        c_2 = float(basis_coefficients_j[b][0])
                        c_3 = float(basis_coefficients_k[c][0])
                        c_4 = float(basis_coefficients_l[d][0])

                        n_1 = ((2 * a_1) / np.pi)**(3/4)
                        n_2 = ((2 * a_2) / np.pi)**(3/4)
                        n_3 = ((2 * a_3) / np.pi)**(3/4)
                        n_4 = ((2 * a_4) / np.pi)**(3/4)

                        x_a = self.basis_set_array[i].get_x()
                        y_a = self.basis_set_array[i].get_y()
                        z_a = self.basis_set_array[i].get_z()
                        x_b = self.basis_set_array[j].get_x()
                        y_b = self.basis_set_array[j].get_y()
                        z_b = self.basis_set_array[j].get_z()
                        x_c = self.basis_set_array[k].get_x()
                        y_c = self.basis_set_array[k].get_y()
                        z_c = self.basis_set_array[k].get_z()
                        x_d = self.basis_set_array[l].get_x()
                        y_d = self.basis_set_array[l].get_y()
                        z_d = self.basis_set_array[l].get_z()
                        x_p = (a_1 * x_a + a_2 * x_b) / (a_1 + a_2)
                        y_p = (a_1 * y_a + a_2 * y_b) / (a_1 + a_2)
                        z_p = (a_1 * z_a + a_2 * z_b) / (a_1 + a_2)
                        x_q = (a_3 * x_c + a_4 * x_d) / (a_3 + a_4)
                        y_q = (a_3 * y_c + a_4 * y_d) / (a_3 + a_4)
                        z_q = (a_3 * z_c + a_4 * z_d) / (a_3 + a_4)

                        r_ab = ((x_a - x_b)**2) + ((y_a - y_b)**2) + ((z_a - z_b)**2)
                        r_cd = ((x_c - x_d)**2) + ((y_c - y_d)**2) + ((z_c - z_d)**2)
                        r_pq = ((x_p - x_q)**2) + ((y_p - y_q)**2) + ((z_p - z_q)**2)

                        s_ab = c_1 * c_2 * n_1 * n_2 * (np.pi / (a_1 + a_2))**(3/2) * np.exp(- a_1 * a_2 * r_ab / (a_1 + a_2))
                        s_cd = c_3 * c_4 * n_3 * n_4 * (np.pi / (a_3 + a_4))**(3/2) * np.exp(- a_3 * a_4 * r_cd / (a_3 + a_4))

                        if r_pq > 0:
                            f_0000 = (1/2) * (np.pi * (a_1 + a_2 + a_3 + a_4) / ((a_1 + a_2) * (a_3 + a_4) * r_pq))**(1/2) * sp.erf(((a_1 + a_2) * (a_3 + a_4) * r_pq / (a_1 + a_2 + a_3 + a_4))**(1/2))
                        else:
                            f_0000 = 1

                        f_mn += 2 * (((a_1 + a_2) * (a_3 + a_4)) / (np.pi * (a_1 + a_2 + a_3 + a_4)))**(1/2) * s_ab * s_cd * f_0000
        return f_mn
