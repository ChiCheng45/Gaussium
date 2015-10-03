import numpy as np


class OverlapIntegral:

    def __init__(self, nuclei_array, file_reader_basis):
        self.nuclei_array = nuclei_array
        self.file_reader_basis = file_reader_basis

    def calculate(self, i, j):
        s_ij = 0
        if i == j:
            return 1.00
        else:
            basis_array_1 = self.file_reader_basis.create_basis_set_array(self.nuclei_array[i].get_name())
            basis_array_2 = self.file_reader_basis.create_basis_set_array(self.nuclei_array[j].get_name())
            x_i = self.nuclei_array[i].get_x()
            y_i = self.nuclei_array[i].get_y()
            z_i = self.nuclei_array[i].get_z()
            x_j = self.nuclei_array[j].get_x()
            y_j = self.nuclei_array[j].get_y()
            z_j = self.nuclei_array[j].get_z()
            for a in range(1, len(basis_array_1)):
                for b in range(1, len(basis_array_2)):
                    if basis_array_1[a][0] == 'S' and basis_array_2[b][0] == 'S':
                        a_1 = float(basis_array_1[a][2])
                        a_2 = float(basis_array_2[b][2])
                        c_1 = float(basis_array_1[a][1])
                        c_2 = float(basis_array_2[b][1])
                        n_1 = ((2 * a_1) / np.pi)**(3/4)
                        n_2 = ((2 * a_2) / np.pi)**(3/4)
                        r_ab = ((x_i - x_j)**2) + ((y_i - y_j)**2) + ((z_i - z_j)**2)
                        s_ij += c_1 * c_2 * n_1 * n_2 * (np.pi / (a_1 + a_2))**(3/2) * np.exp(-a_1 * a_2 * r_ab / (a_1 + a_2))
                    else:
                        s_ij += 0
            return s_ij
