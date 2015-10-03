from src.main import Coulomb, CoulombTotal, FileInputNuclei, FileInputBasis, Matrix, OverlapIntegral, \
    KineticEnergyIntegral, NuclearAttractionIntegral, TwoElectronRepulsion, DensityMatrix, SCFProcedure, TotalEnergy
import numpy as np
import time

if __name__ == '__main__':
    start = time.clock()
    print('*********************************************************************************************************')
    print('\nA BASIC QUANTUM CHEMICAL PROGRAM IN PYTHON\n\n')

    file_reader_nuclei = FileInputNuclei('HeH+.mol')
    file_reader_basis = FileInputBasis('STO-3G.gbs')
    nuclei_array = file_reader_nuclei.create_nuclei_array()
    nuclei_name_array = []
    for a in range(0, len(nuclei_array)):
        nuclei_name_array.append(nuclei_array[a].get_name())
    print(nuclei_name_array)

    coulomb_total = CoulombTotal(Coulomb, nuclei_array)
    coulomb_law_matrix = coulomb_total.calculate_total_electric_potential_energy()
    nuclear_repulsion_energy = coulomb_law_matrix.sum() / 2
    print('\nCoulombs Law Matrix')
    print(coulomb_law_matrix)
    print('Total Nuclear-Nuclear Potential Energy: ' + str(nuclear_repulsion_energy) + ' a.u.')

    matrix = Matrix(nuclei_array)
    s_matrix = matrix.create_matrix(OverlapIntegral(nuclei_array, file_reader_basis))
    print('\nOrbital Overlap Matrix')
    print(s_matrix)
    t_matrix = matrix.create_matrix(KineticEnergyIntegral(nuclei_array, file_reader_basis))
    print('\nKinetic Energy Matrix')
    print(t_matrix)
    v_matrix = matrix.create_matrix(NuclearAttractionIntegral(nuclei_array, file_reader_basis))
    print('\nNuclear Potential Energy Matrix')
    print(v_matrix)
    tei_matrix = matrix.create_matrix(TwoElectronRepulsion(nuclei_array, file_reader_basis))
    print('\nTwo Electron Repulsion Energy Matrix')
    print('------------[11, 12, 22]------------')
    print(tei_matrix.item(1, 0))

    print('\n*********************************************************************************************************')
    print('\nH_CORE AND TRANSFORMATION MATRIX\n')

    h_core_matrix = t_matrix + v_matrix
    print('\nH_CORE')
    print(h_core_matrix)

    s_matrix_eigenvalues = np.linalg.eig(s_matrix)[0]
    s_matrix_unitary = np.linalg.eig(s_matrix)[1]
    x_canonical = s_matrix_unitary * np.diag(s_matrix_eigenvalues ** (-1/2))
    print('\nTRANSFORMATION MATRIX')
    print(x_canonical)

    print('\n*********************************************************************************************************')
    print('\nDENSITY MATRIX INITIAL GUESS\n')

    """Create a initial guess for the density matrix by turning off all two-electron interaction and solving for the
    orbital coefficients. The two-electron parts are then turned back on during the SCF procedure.
    """

    orthonormal_h_matrix = x_canonical.T * h_core_matrix * x_canonical
    orbital_energy_matrix = np.diag(np.linalg.eig(orthonormal_h_matrix)[0])
    print('\nORBITAL ENERGY EIGENVALUES')
    print(orbital_energy_matrix)

    """The negative signs that numpy give for the orbital coefficients seems strange but it is because the sign of the
    orbitals and its wave-function can either have positive or negative values, the square of the wave-function, the
    electron density, will always end up being positive valued. The elements in the density matrix will made by
    multiplying two of the orbital coefficients together anyway
    """

    orbital_coefficients_orthonormal = np.linalg.eig(orthonormal_h_matrix)[1]
    orbital_coefficients = x_canonical * orbital_coefficients_orthonormal
    print('\nORBITAL COEFFICIENTS')
    print(orbital_coefficients)

    density_matrix = DensityMatrix(orbital_coefficients)
    p_matrix = matrix.create_matrix(density_matrix)
    print('\nDENSITY MATRIX')
    print(p_matrix)

    print('\n*********************************************************************************************************')
    print('\nBEGIN SCF PROCEDURE')
    total_energy = TotalEnergy()
    scf_procedure = SCFProcedure(tei_matrix, h_core_matrix, x_canonical, matrix, total_energy, nuclear_repulsion_energy)
    scf_procedure.begin_scf(p_matrix)

    print('\n*********************************************************************************************************')
    print('\nTime Taken: ' + str(time.clock() - start) + 's')
    print("\nWhat I cannot create I cannot understand - Richard Feynman\n")
