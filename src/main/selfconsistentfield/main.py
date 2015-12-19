from src.main.common import Matrix
from src.main.fileinput import FileInputBasis, FileInputNuclei
from src.main.matrixelements import KineticEnergyElement, NuclearAttractionElement, OrbitalOverlapElement, TwoElectronRepulsionElementCook, TwoElectronRepulsionElementOS, TwoElectronRepulsionElementHGP
from src.main.selfconsistentfield import SCFProcedure, CoulombsLawArray
import numpy as np
import time


if __name__ == '__main__':
    np.set_printoptions(linewidth=100000, threshold=np.inf)

    start = time.clock()
    print('***********************************************************************************************************')
    print('\nA BASIC QUANTUM CHEMICAL PROGRAM IN PYTHON\n\n')

    nuclei_array, electrons = FileInputNuclei('C2H4.mol').create_nuclei_array_and_electron_count()
    basis_set_array = FileInputBasis('3-21G.gbs', nuclei_array).create_basis_set_array()

    nuclei_name_list = [x.element for x in nuclei_array]
    print(nuclei_name_list)

    coulomb_law_array = CoulombsLawArray.calculate_total_electric_potential_energy(nuclei_array)
    nuclear_repulsion = coulomb_law_array.sum() / 2
    print('\nNUCLEAR REPULSION ARRAY')
    print(coulomb_law_array)

    print('\n*********************************************************************************************************')
    print('\nMATRICES\n')

    matrix = Matrix(len(basis_set_array))

    print('\nORBITAL OVERLAP MATRIX')
    s_matrix = matrix.create_matrix(OrbitalOverlapElement(basis_set_array))
    print(s_matrix)

    print('\nKINETIC ENERGY MATRIX')
    t_matrix = matrix.create_matrix(KineticEnergyElement(basis_set_array))
    print(t_matrix)

    print('\nNUCLEAR POTENTIAL ENERGY MATRIX')
    v_matrix = matrix.create_matrix(NuclearAttractionElement(nuclei_array, basis_set_array))
    print(v_matrix)

    print('\nCORE HAMILTONIAN MATRIX')
    h_core_matrix = t_matrix + v_matrix
    print(h_core_matrix)

    print('\nTRANSFORMATION MATRIX')
    s_matrix_eigenvalues, s_matrix_unitary = np.linalg.eigh(s_matrix)
    sort = np.argsort(s_matrix_eigenvalues)
    s_matrix_eigenvalues = np.array(s_matrix_eigenvalues)[sort]
    s_matrix_unitary = s_matrix_unitary[:, sort]
    s_matrix_eigenvalues = [x**(-1/2) for x in s_matrix_eigenvalues]
    x_canonical = s_matrix_unitary * np.diag(s_matrix_eigenvalues)
    print(x_canonical)

    print('\n*********************************************************************************************************')
    print('\nINITIAL GUESS\n')
    """
    Create a initial guess for the density matrix from the solutions of the orbital coefficients with all two-electron
    interactions turned off. The two-electron parts are then turned back on during the SCF procedure.
    """

    print('\nORBITAL ENERGY EIGENVALUES')
    orthonormal_h_matrix = np.transpose(x_canonical) * h_core_matrix * x_canonical
    eigenvalues, eigenvectors = np.linalg.eigh(orthonormal_h_matrix)
    sort = np.argsort(eigenvalues)
    eigenvalues = np.array(eigenvalues)[sort]
    print(eigenvalues)

    print('\nORBITAL COEFFICIENTS')
    eigenvectors = eigenvectors[:, sort]
    orbital_coefficients = x_canonical * eigenvectors
    print(orbital_coefficients)

    print('\n*********************************************************************************************************')
    print('\nBEGIN SCF PROCEDURE')
    """
    There are two methods to produce the repulsion dictionary. One uses multiprocessing and gives some speed up on
    larger basis sets molecule. For smaller basis set molecules its better to use the normal single core method. My
    computer is a dual core with Hyper-threading and I have found that running four processes gives the most performance
    boost for the molecules I have tested. Doing this actually maxes my cpu out during this part of the calculation.
    """

    # repulsion_dictionary = TwoElectronRepulsionElementCook(basis_set_array).store_parallel(4)
    repulsion_dictionary = TwoElectronRepulsionElementOS(basis_set_array).store_parallel(4)
    # repulsion_dictionary = TwoElectronRepulsionElementHGP(basis_set_array).store_parallel(4)

    scf_procedure = SCFProcedure(h_core_matrix, x_canonical, matrix, electrons, repulsion_dictionary)
    electron_energy = scf_procedure.begin_scf(orbital_coefficients)

    print('TOTAL NUCLEAR REPULSION ENERGY: ' + str(nuclear_repulsion) + ' a.u.')
    print('TOTAL ENERGY: ' + str(electron_energy + nuclear_repulsion) + ' a.u.')

    print('\n*********************************************************************************************************')
    print('\nTime Taken: ' + str(time.clock() - start) + 's')
    print("\nWhat I cannot create I cannot understand - Richard Feynman\n")
