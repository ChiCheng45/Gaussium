from src.main.objects import expand_basis_set
from src.main.objects import Nuclei
import os, sys, re


def read_basis_set_file(file_input_basis, nuclei_array):
    file_input_basis = os.path.join(sys.path[1], 'basissets\\' + file_input_basis)
    basis_array = []
    for a in range(len(nuclei_array)):
        file = open(file_input_basis, 'r')
        lines = file.read().replace('\n', ':')
        file.close()

        regex = nuclei_array[a].element + '.*?#'
        lines = ' '.join(lines.split())
        lines += '#'
        lines = re.search(regex, lines).group(0)
        lines = re.split(':', lines.replace(': ', ':'))
        lines = lines[1:len(lines) - 2]

        i = 0
        input1 = input2 = []
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
        basis_array_from_fact = expand_basis_set(input1, nuclei_array[a].coordinates)
        basis_array += basis_array_from_fact

    return basis_array


def read_mol_file(file_input_mol):
    file_input_mol = os.path.join(sys.path[1], 'molfiles\\' + file_input_mol)
    nuclei_array = []
    total_nuclei_charge = 0
    with open(file_input_mol, 'r') as file:
        lines = file.readlines()
        molecular_charge = int(lines[0].split()[0])
        multiplicity = int(lines[0].split()[1])
        for a in range(1, len(lines)):
            line = lines[a]
            array = line.split()
            total_nuclei_charge += int(array[1])
            coordinates = (float(array[3]), float(array[4]), float(array[5]))
            nuclei = Nuclei(array[0], float(array[1]), float(array[2]), coordinates)
            nuclei_array.append(nuclei)
    electron_count = total_nuclei_charge - molecular_charge
    return nuclei_array, electron_count, multiplicity
