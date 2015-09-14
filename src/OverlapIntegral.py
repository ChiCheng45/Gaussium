class OverlapIntegral:

    def calculate_integral(self, nuclei_array, file_reader_basis, i, j):
        s_ij = 0
        if i == j:
            return 1.00
        else:
            basis_array_1 = file_reader_basis.create_basis_set_array(nuclei_array[i].get_name())
            basis_array_2 = file_reader_basis.create_basis_set_array(nuclei_array[j].get_name())
            for a in range(1, len(basis_array_1)):
                for b in range(1, len(basis_array_2)):
                    if basis_array_1[a][0] == 'S' and basis_array_2[b][0] == 'S':
                        s_ij += float(basis_array_1[a][2]) * float(basis_array_2[b][2])
                    else:
                        s_ij += 0
            return s_ij