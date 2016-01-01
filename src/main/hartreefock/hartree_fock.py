from src.main.common import Matrix
from src.main.matrixelements import KineticEnergyElement, NuclearAttractionElement, OrbitalOverlapElement, TwoElectronRepulsionElementCook, TwoElectronRepulsionElementOS, TwoElectronRepulsionElementHGP
from src.main.hartreefock import SCF, LinearAlgebra
import time


class HartreeFock:

    @staticmethod
    def begin(nuclei_array, electrons, multiplicity, basis_set_array):
        print('\n*****************************************************************************************************')
        create = Matrix(len(basis_set_array)).create_matrix
        orbital_overlap = create(OrbitalOverlapElement(basis_set_array))
        kinetic_energy = create(KineticEnergyElement(basis_set_array))
        nuclear_potential = create(NuclearAttractionElement(nuclei_array, basis_set_array))
        core_hamiltonian = kinetic_energy + nuclear_potential
        transformation_matrix = LinearAlgebra.create_transformation_matrix(orbital_overlap)

        print('\nMATRICES\n')
        print('\nORBITAL OVERLAP MATRIX')
        print(orbital_overlap)
        print('\nKINETIC ENERGY MATRIX')
        print(kinetic_energy)
        print('\nNUCLEAR POTENTIAL ENERGY MATRIX')
        print(nuclear_potential)
        print('\nCORE HAMILTONIAN MATRIX')
        print(core_hamiltonian)
        print('\nTRANSFORMATION MATRIX')
        print(transformation_matrix)

        print('\n*****************************************************************************************************')
        diagonalize = LinearAlgebra(transformation_matrix).diagonalize
        orbital_energies, orbital_coefficients = diagonalize(core_hamiltonian)

        print('\nINITIAL GUESS\n')
        print('\nORBITAL ENERGY EIGENVALUES')
        print(orbital_energies)
        print('\nORBITAL COEFFICIENTS')
        print(orbital_coefficients)

        print('\n*****************************************************************************************************')
        print('\nBEGIN TWO ELECTRON REPULSION CALCULATION')
        start_repulsion = time.clock()
        # repulsion_dictionary = TwoElectronRepulsionElementCook(basis_set_array).store_parallel(4)
        # repulsion_dictionary = TwoElectronRepulsionElementHGP(basis_set_array).store_parallel(4)
        repulsion_dictionary = TwoElectronRepulsionElementOS(basis_set_array).store_parallel(4)
        print('TIME TAKEN: ' + str(time.clock() - start_repulsion) + 's\n')

        scf_procedure = SCF(core_hamiltonian, diagonalize, create, electrons, repulsion_dictionary)
        return scf_procedure, orbital_coefficients, multiplicity, repulsion_dictionary

    @classmethod
    def restricted(cls, nuclei_array, electrons, multiplicity, basis_set_array):
        scf_procedure, orbital_coefficients, multiplicity, repulsion = cls.begin(nuclei_array, electrons, multiplicity, basis_set_array)

        print('\nBEGIN SCF PROCEDURE')
        start = time.clock()
        electron_energy, eigenvalues, orbital_coefficients = scf_procedure.begin_rhf(orbital_coefficients)
        print('TIME TAKEN: ' + str(time.clock() - start) + 's\n')

        print('\nORBITAL ENERGY EIGENVALUES')
        print(eigenvalues)
        print('\nORBITAL COEFFICIENTS')
        print(orbital_coefficients, end='\n\n')
        return electron_energy, eigenvalues, orbital_coefficients, repulsion

    @classmethod
    def unrestricted(cls, nuclei_array, electrons, multiplicity, basis_set_array):
        scf_procedure, orbital_coefficients, multiplicity, repulsion = cls.begin(nuclei_array, electrons, multiplicity, basis_set_array)

        print('\nBEGIN SCF PROCEDURE')
        start = time.clock()
        electron_energy, eigenvalues_alpha, eigenvalues_beta, orbital_coefficients_alpha, orbital_coefficients_beta = scf_procedure.begin_uhf(orbital_coefficients, multiplicity)
        print('TIME TAKEN: ' + str(time.clock() - start) + 's\n')

        print('\nALPHA ORBITAL ENERGY EIGENVALUES')
        print(eigenvalues_alpha)
        print('\nBETA ORBITAL ENERGY EIGENVALUES')
        print(eigenvalues_beta)
        print('\nALPHA ORBITAL COEFFICIENTS')
        print(orbital_coefficients_alpha, end='\n\n')
        print('\nBETA ORBITAL COEFFICIENTS')
        print(orbital_coefficients_beta, end='\n\n')
        return electron_energy, eigenvalues_alpha, eigenvalues_beta, orbital_coefficients_alpha, orbital_coefficients_beta, repulsion