import os
from src.CoulombsLaw import Coulomb
from src.FileInputNuclei import FileInputNuclei
from src.FileInputBasis import FileInputBasis
from src.CreateMatrix import Matrix
from src.OverlapIntegral import OverlapIntegral


def main():
    print("\nWhat I cannot create I cannot understand - Richard Feynman\n")

    file_reader_nuclei = FileInputNuclei('HeH+.mol')
    file_reader_basis = FileInputBasis('STO-3G-edited.gbs')

    nuclei_array = file_reader_nuclei.create_nuclei_array()
    nuclei_name_array = []
    for a in range(0, len(nuclei_array)):
        nuclei_name_array.append(nuclei_array[a].get_name())

    coulomb = Coulomb(nuclei_array)
    print('Nuclear-Nuclear Potential Energy: ' + str(coulomb.calculate_total_electric_potential_energy()) + ' E_h \n')

    matrix = Matrix(nuclei_array, file_reader_basis)
    s_matrix = matrix.create_matrix(OverlapIntegral())
    print(nuclei_name_array)
    print(s_matrix)

main()
