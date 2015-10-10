import os
import sys
from src.main.basis import Basis


class FileInputBasis:

    def __init__(self, file_input_basis, nuclei_array):
        self.file_input_mol = os.path.join(sys.path[1], 'basisSetFiles\\' + file_input_basis)
        self.nuclei_array = nuclei_array

    def create_basis_set_array(self):
        basis_array = []
        for a in range(0, len(self.nuclei_array)):
            i = j = 0
            nuclei = self.nuclei_array[a]
            coefficients_array = []
            with open(self.file_input_mol, 'r') as file:
                lines = file.readlines()
                for b in range(0, len(lines)):
                    line = lines[b]
                    if nuclei.get_name() in line:
                        i = 1
                    if i == 1:
                        if '#' in line:
                            if nuclei.get_name() not in line:
                                break
                        else:
                            if any(letter in line for letter in ('S', 'L', 'P', 'D')) or line == '\n':
                                if j == 1:
                                    basis = Basis(nuclei.get_name(), nuclei.get_y(), nuclei.get_x(), nuclei.get_z(), function_type, coefficients_array)
                                    basis_array.append(basis)
                                    j = 0
                                if line != '\n':
                                    coefficients_array = []
                                    function_type = line.split()[0]
                            else:
                                j = 1
                                float_array = [float(x) for x in line.split()]
                                coefficients_array.append(float_array)
                                if b + 1 == len(lines):
                                    basis = Basis(nuclei.get_name(), nuclei.get_y(), nuclei.get_x(), nuclei.get_z(), function_type, coefficients_array)
                                    basis_array.append(basis)
            file.close()
        return basis_array
