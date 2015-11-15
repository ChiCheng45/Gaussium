from src.main.common import Matrix
from src.main.fileinput import FileInputBasis, FileInputNuclei
from src.main.matrixelements import KineticEnergyElement, NuclearAttractionElement, OrbitalOverlapElement, TwoElectronRepulsionElement
from src.main.selfconsistentfield import SCFProcedure, CoulombsLawArray
import numpy as np
import time


if __name__ == '__main__':
    np.set_printoptions(linewidth=10000)

    start = time.clock()
    print('*********************************************************************************************************')
    print('\nA BASIC QUANTUM CHEMICAL PROGRAM IN PYTHON\n\n')

    nuclei_array, electrons = FileInputNuclei('HeH+.mol').create_nuclei_array_and_electron_count()
    basis_set_array = FileInputBasis('6-311+GPP.gbs', nuclei_array).create_basis_set_array()

    nuclei_name_list = [x.element for x in nuclei_array]
    print(nuclei_name_list)

    coulomb_law_array = CoulombsLawArray.calculate_total_electric_potential_energy(nuclei_array)
    nuclear_repulsion = coulomb_law_array.sum() / 2
    print('\nNUCLEAR REPULSION ARRAY')
    print(coulomb_law_array)

    print('\n*********************************************************************************************************')
    print('\nMATRICES\n')

    matrix = Matrix(len(basis_set_array))

    s_matrix = matrix.create_matrix(OrbitalOverlapElement(basis_set_array))
    print('\nORBITAL OVERLAP MATRIX')
    print(s_matrix)

    t_matrix = matrix.create_matrix(KineticEnergyElement(basis_set_array))
    print('\nKINETIC ENERGY MATRIX')
    print(t_matrix)

    v_matrix = matrix.create_matrix(NuclearAttractionElement(nuclei_array, basis_set_array))
    print('\nNUCLEAR POTENTIAL ENERGY MATRIX')
    print(v_matrix)

    h_core_matrix = t_matrix + v_matrix
    print('\nCORE HAMILTONIAN MATRIX')
    print(h_core_matrix)

    s_matrix_eigenvalues, s_matrix_unitary = np.linalg.eig(s_matrix)
    sort = s_matrix_eigenvalues.argsort()[::1]
    s_matrix_eigenvalues = s_matrix_eigenvalues[sort]
    s_matrix_eigenvalues = [x**(-1/2) for x in s_matrix_eigenvalues]
    s_matrix_unitary = s_matrix_unitary[:, sort]
    x_canonical = s_matrix_unitary * np.diag(s_matrix_eigenvalues)
    print('\nTRANSFORMATION MATRIX')
    print(x_canonical)

    print('\n*********************************************************************************************************')
    print('\nINITIAL GUESS\n')

    """
    Create a initial guess for the density matrix from the solutions of the orbital coefficients with all two-electron
    interactions turned off. The two-electron parts are then turned back on during the SCF procedure.
    """

    orthonormal_h_matrix = x_canonical.T * h_core_matrix * x_canonical
    eigenvalues, eigenvectors = np.linalg.eig(orthonormal_h_matrix)
    sort = eigenvalues.argsort()[::1]
    eigenvalues = eigenvalues[sort]
    eigenvectors = eigenvectors[:, sort]

    orbital_energy_matrix = np.diag(eigenvalues)
    print('\nORBITAL ENERGY EIGENVALUES')
    print(orbital_energy_matrix)

    orbital_coefficients = x_canonical * eigenvectors
    print('\nORBITAL COEFFICIENTS')
    print(orbital_coefficients)

    print('\n*********************************************************************************************************')
    print('\nBEGIN SCF PROCEDURE')

    repulsion_dictionary = TwoElectronRepulsionElement(basis_set_array).store_integrals()

    scf_procedure = SCFProcedure(h_core_matrix, x_canonical, matrix, electrons)
    electron_energy = scf_procedure.begin_scf(orbital_coefficients, repulsion_dictionary)

    print('TOTAL NUCLEAR REPULSION ENERGY: ' + str(nuclear_repulsion) + ' a.u.')
    print('TOTAL ENERGY: ' + str(electron_energy + nuclear_repulsion) + ' a.u.')

    print('\n*********************************************************************************************************')
    print('\nTime Taken: ' + str(time.clock() - start) + 's')
    print("\nWhat I cannot create I cannot understand - Richard Feynman\n")
