import numpy as np
import scipy.special as sp

# some really bad code below, sorry

class TwoElectronRepulsion:

    def __init__(self, nuclei_array, file_reader_basis):
        self.nuclei_array = nuclei_array
        self.file_reader_basis = file_reader_basis

    def calculate(self, i, j):
        basis_array_1 = self.file_reader_basis.create_basis_set_array(self.nuclei_array[i].get_name())
        basis_array_2 = self.file_reader_basis.create_basis_set_array(self.nuclei_array[j].get_name())

        array = [basis_array_2, basis_array_1]
        nuclei_array_2 = [self.nuclei_array[j], self.nuclei_array[i]]
        matrix = []
        for a in range(0, 2):
            for b in range(0, 2):
                if not (a == 0 and b == 1):
                    row = []
                    for c in range(0, 2):
                        for d in range(0, 2):
                            if not (c == 0 and d == 1):
                                f_mn = 0
                                for e in range(1, len(array[a])):
                                    for f in range(1, len(array[b])):
                                        for g in range(1, len(array[c])):
                                            for h in range(1, len(array[d])):

                                                a_1 = float(array[a][e][2])
                                                a_2 = float(array[b][f][2])
                                                a_3 = float(array[c][g][2])
                                                a_4 = float(array[d][h][2])

                                                c_1 = float(array[a][e][1])
                                                c_2 = float(array[b][f][1])
                                                c_3 = float(array[c][g][1])
                                                c_4 = float(array[d][h][1])

                                                n_1 = ((2 * a_1) / np.pi)**(3/4)
                                                n_2 = ((2 * a_2) / np.pi)**(3/4)
                                                n_3 = ((2 * a_3) / np.pi)**(3/4)
                                                n_4 = ((2 * a_4) / np.pi)**(3/4)

                                                x_a = nuclei_array_2[a].get_x()
                                                y_a = nuclei_array_2[a].get_y()
                                                z_a = nuclei_array_2[a].get_z()
                                                x_b = nuclei_array_2[b].get_x()
                                                y_b = nuclei_array_2[b].get_y()
                                                z_b = nuclei_array_2[b].get_z()
                                                x_c = nuclei_array_2[c].get_x()
                                                y_c = nuclei_array_2[c].get_y()
                                                z_c = nuclei_array_2[c].get_z()
                                                x_d = nuclei_array_2[d].get_x()
                                                y_d = nuclei_array_2[d].get_y()
                                                z_d = nuclei_array_2[d].get_z()
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
                                row.append(f_mn)
                    matrix.append(row)
        return np.matrix(matrix)
