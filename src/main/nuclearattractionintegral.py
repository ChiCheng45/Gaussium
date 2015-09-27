import numpy as np
import scipy.special as sp


class NuclearAttractionIntegral:

    def calculate_integral(self, nuclei_array, file_reader_basis, i, j):
        v_ij = 0
        basis_array_1 = file_reader_basis.create_basis_set_array(nuclei_array[i].get_name())
        basis_array_2 = file_reader_basis.create_basis_set_array(nuclei_array[j].get_name())
        x_i = nuclei_array[i].get_x()
        y_i = nuclei_array[i].get_y()
        z_i = nuclei_array[i].get_z()
        x_j = nuclei_array[j].get_x()
        y_j = nuclei_array[j].get_y()
        z_j = nuclei_array[j].get_z()
        for a in range(1, len(basis_array_1)):
            for b in range(1, len(basis_array_2)):
                if basis_array_1[a][0] == 'S' and basis_array_2[b][0] == 'S':
                    a_1 = float(basis_array_1[a][2])
                    a_2 = float(basis_array_2[b][2])
                    d_1 = float(basis_array_1[a][1])
                    d_2 = float(basis_array_2[b][1])
                    n_1 = ((2 * a_1) / np.pi)**(3/4)
                    n_2 = ((2 * a_2) / np.pi)**(3/4)
                    r_ij = ((x_i - x_j)**2) + ((y_i - y_j)**2) + ((z_i - z_j)**2)
                    s_ij = d_1 * d_2 * n_1 * n_2 * (np.pi / (a_1 + a_2))**(3/2) * np.exp(- a_1 * a_2 * r_ij / (a_1 + a_2))
                    for k in range(0, len(nuclei_array)):
                        x_k = nuclei_array[k].get_x()
                        y_k = nuclei_array[k].get_y()
                        z_k = nuclei_array[k].get_z()
                        c_k = nuclei_array[k].get_charge()
                        x_3 = (a_1 * x_i + a_2 * x_j) / (a_1 + a_2)
                        y_3 = (a_1 * y_i + a_2 * y_j) / (a_1 + a_2)
                        z_3 = (a_1 * z_i + a_2 * z_j) / (a_1 + a_2)
                        r_3k = ((x_k - x_3)**2) + ((y_k - y_3)**2) + ((z_k - z_3)**2)
                        if r_3k > 0:
                            f_ij = (1/2) * (np.pi / ((a_1 + a_2) * r_3k))**(1/2) * sp.erf(((a_1 + a_2) * r_3k)**(1/2))
                        else:
                            f_ij = 1
                        v_ij += - 2 * c_k * np.sqrt((a_1 + a_2) / np.pi) * s_ij * f_ij
                else:
                    v_ij += 0
        return v_ij

