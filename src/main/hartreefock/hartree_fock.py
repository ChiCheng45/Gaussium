from src.main.common import Matrix, CoulombsLawArray
from src.main.fileinput import FileInputBasis, FileInputNuclei
from src.main.matrixelements import KineticEnergyElement, NuclearAttractionElement, OrbitalOverlapElement, TwoElectronRepulsionElementCook, TwoElectronRepulsionElementOS, TwoElectronRepulsionElementHGP
from src.main.hartreefock import SCF
import numpy as np
import time


class HartreeFock:

    @staticmethod
    def begin(mol, basis, method):
        np.set_printoptions(linewidth=100000, threshold=np.inf)

        start = time.clock()
        print('*******************************************************************************************************')
        print('\nA BASIC QUANTUM CHEMICAL PROGRAM IN PYTHON\n\n')

        nuclei_array, electrons, multiplicity = FileInputNuclei(mol).create_nuclei_array_and_electron_count()
        basis_set_array = FileInputBasis(basis, nuclei_array).create_basis_set_array()

        nuclei_name_list = [x.element for x in nuclei_array]
        print(nuclei_name_list)

        coulomb_law_array = CoulombsLawArray.calculate_total_electric_potential_energy(nuclei_array)
        nuclear_repulsion = coulomb_law_array.sum() / 2
        print('\nNUCLEAR REPULSION ARRAY')
        print(coulomb_law_array)

        print('\n*****************************************************************************************************')
        print('\nMATRICES\n')

        create = Matrix(len(basis_set_array)).create_matrix

        print('\nORBITAL OVERLAP MATRIX')
        orbital_overlap = create(OrbitalOverlapElement(basis_set_array))
        print(orbital_overlap)

        print('\nKINETIC ENERGY MATRIX')
        kinetic_energy = create(KineticEnergyElement(basis_set_array))
        print(kinetic_energy)

        print('\nNUCLEAR POTENTIAL ENERGY MATRIX')
        nuclear_potential = create(NuclearAttractionElement(nuclei_array, basis_set_array))
        print(nuclear_potential)

        print('\nCORE HAMILTONIAN MATRIX')
        core_hamiltonian = kinetic_energy + nuclear_potential
        print(core_hamiltonian)

        print('\nTRANSFORMATION MATRIX')
        s_matrix_eigenvalues, s_matrix_unitary = np.linalg.eigh(orbital_overlap)
        sort = np.argsort(s_matrix_eigenvalues)
        s_matrix_eigenvalues = np.array(s_matrix_eigenvalues)[sort]
        s_matrix_unitary = s_matrix_unitary[:, sort]
        s_matrix_eigenvalues = [x**(-1/2) for x in s_matrix_eigenvalues]
        x_canonical = s_matrix_unitary * np.diag(s_matrix_eigenvalues)
        print(x_canonical)

        print('\n*****************************************************************************************************')
        print('\nINITIAL GUESS\n')

        print('\nORBITAL ENERGY EIGENVALUES')
        orthonormal_h_matrix = np.transpose(x_canonical) * core_hamiltonian * x_canonical
        eigenvalues, eigenvectors = np.linalg.eigh(orthonormal_h_matrix)
        sort = np.argsort(eigenvalues)
        eigenvalues = np.array(eigenvalues)[sort]
        eigenvectors = eigenvectors[:, sort]
        orbital_coefficients = x_canonical * eigenvectors
        print(eigenvalues)

        print('\nORBITAL COEFFICIENTS')
        print(orbital_coefficients)

        print('\n*****************************************************************************************************')

        print('\nBEGIN TWO ELECTRON REPULSION CALCULATION')
        start_repulsion = time.clock()
        # repulsion_dictionary = TwoElectronRepulsionElementCook(basis_set_array).store_parallel(4)
        # repulsion_dictionary = TwoElectronRepulsionElementHGP(basis_set_array).store_parallel(4)
        repulsion_dictionary = TwoElectronRepulsionElementOS(basis_set_array).store_parallel(4)
        print('TIME TAKEN: ' + str(time.clock() - start_repulsion) + 's\n')

        print('\nBEGIN SCF PROCEDURE')
        scf_procedure = SCF(core_hamiltonian, x_canonical, create, electrons, repulsion_dictionary)
        if method == 'RHF':
            electron_energy = scf_procedure.begin_rhf(orbital_coefficients)
        elif method == 'UHF':
            electron_energy = scf_procedure.begin_uhf(orbital_coefficients, multiplicity)
        else:
            electron_energy = 0
        print('TOTAL NUCLEAR REPULSION ENERGY: ' + str(nuclear_repulsion) + ' a.u.')
        print('TOTAL ENERGY: ' + str(electron_energy + nuclear_repulsion) + ' a.u.')

        print('\n*****************************************************************************************************')
        print('\nTIME TAKEN: ' + str(time.clock() - start) + 's')
        print("\nWhat I cannot create I cannot understand - Richard Feynman\n")
