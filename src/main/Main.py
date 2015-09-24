import time
from src.main import Coulomb, CoulombTotal, FileInputNuclei, FileInputBasis, Matrix, OverlapIntegral, KineticEnergyIntegral


def main():
    print("\nWhat I cannot create I cannot understand - Richard Feynman\n")

    start = time.clock()
    file_reader_nuclei = FileInputNuclei('HeH+.mol')
    file_reader_basis = FileInputBasis('STO-3G-edited.gbs')

    nuclei_array = file_reader_nuclei.create_nuclei_array()
    nuclei_name_array = []
    for a in range(0, len(nuclei_array)):
        nuclei_name_array.append(nuclei_array[a].get_name())

    print(nuclei_name_array)

    coulomb_total = CoulombTotal(Coulomb, nuclei_array)
    coulomb_law_matrix = coulomb_total.calculate_total_electric_potential_energy()
    print('\nCoulombs Law Matrix')
    print(coulomb_law_matrix)
    print('Total Nuclear Potential Energy: ' + str(coulomb_law_matrix.sum() / 2) + ' E_h')

    matrix = Matrix(nuclei_array, file_reader_basis)
    s_matrix = matrix.create_matrix(OverlapIntegral())
    print('\nOrbital Overlap Matrix')
    print(s_matrix)

    t_matrix = matrix.create_matrix(KineticEnergyIntegral())
    print('\nKinetic Energy Matrix')
    print(t_matrix)
    print('\nTime Taken: ' + str(time.clock() - start) + 's')

main()
