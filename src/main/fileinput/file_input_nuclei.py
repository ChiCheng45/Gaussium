import os
import sys

from src.main.objects import Nuclei


class FileInputNuclei:

    def __init__(self, file_input_mol):
        self.file_input_mol = os.path.join(sys.path[1], 'molfiles\\' + file_input_mol)

    def create_nuclei_array(self):
        nuclei_array = []
        file = open(self.file_input_mol, 'r')
        for line in file:
            array = line.split()
            nuclei = Nuclei(array)
            nuclei_array.append(nuclei)
        file.close()
        return nuclei_array
