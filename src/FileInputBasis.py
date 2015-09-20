import os


class FileInputBasis:

    def __init__(self, file_input_basis):
        self.file_input_basis = os.path.abspath('./basisSetFiles/' + file_input_basis)

    def create_basis_set_array(self, nuclei_name):
        i = j = 0
        basis_array = [nuclei_name]
        with open(self.file_input_basis, 'r') as file:
            lines = file.readlines()
            for a in range(0, len(lines)):
                line = lines[a]
                if nuclei_name in line:
                    i = 1
                if i == 1:
                    if '#' in line and nuclei_name not in line:
                        break
                    if '#' not in line and any(letter in line for letter in ('S', 'L', 'D')):
                        j = 1
                        function_type = line.split()[0]
                    if j == 1 and not any(letter in line for letter in ('S', 'L', 'D')) and line != '\n':
                        function_array = [function_type]
                        for b in range(0, len(line.split())):
                            function_array.append(line.split()[b])
                        basis_array.append(function_array)
        file.close()
        return basis_array
