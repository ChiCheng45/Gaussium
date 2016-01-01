import os
import sys
import re
from src.main.objects import PrimitiveBasisFactory


class FileInputBasis:

    @staticmethod
    def create_basis_set_array(file_input_basis, nuclei_array):
        file_input_basis = os.path.join(sys.path[1], 'basissets\\' + file_input_basis)
        basis_array = []
        for a in range(len(nuclei_array)):
            nuclei = nuclei_array[a]
            regex = nuclei.element + '.*?#'
            file = open(file_input_basis, 'r')
            lines = file.read().replace('\n', ':')
            file.close()
            lines = ' '.join(lines.split())
            lines = re.search(regex, lines).group(0)
            lines = re.split(':', lines.replace(': ', ':'))
            lines = lines[1:len(lines) - 2]
            i = 0
            input1 = []
            for line in lines:
                if any(letter in line for letter in ('S', 'L', 'P', 'D')):
                    if i == 1:
                        input1.append(input2)
                        input2 = [line[0]]
                    else:
                        i = 1
                        input2 = [line[0]]
                else:
                    input2.append([float(b) for b in line.split()])
            input1.append(input2)
            basis_array_from_fact = PrimitiveBasisFactory.expand_basis(input1, nuclei.coordinates)
            basis_array += basis_array_from_fact
        return basis_array
