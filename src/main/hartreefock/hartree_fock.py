from src.main.common import Matrix
from src.main.matrixelements import KineticEnergyElement, NuclearAttractionElement, OrbitalOverlapElement, TwoElectronRepulsionElementCook, TwoElectronRepulsionElementOS, TwoElectronRepulsionElementHGP
from src.main.hartreefock import SelfConsistentField, LinearAlgebra
import time


class HartreeFock:

    def __init__(self, nuclei_array, electrons, multiplicity, basis_set_array):
        self.nuclei_array = nuclei_array
        self.electrons = electrons
        self.multiplicity = multiplicity
        self.basis_set_array = basis_set_array

    def begin(self):
        print('\n*****************************************************************************************************')
        create_matrix = Matrix(len(self.basis_set_array)).create_matrix
        orbital_overlap = create_matrix(OrbitalOverlapElement(self.basis_set_array))
        kinetic_energy = create_matrix(KineticEnergyElement(self.basis_set_array))
        nuclear_potential = create_matrix(NuclearAttractionElement(self.nuclei_array, self.basis_set_array))
        core_hamiltonian = kinetic_energy + nuclear_potential
        transformation_matrix = LinearAlgebra.transformation_matrix(orbital_overlap)

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
        print('\nINITIAL GUESS\n')
        diagonalize = LinearAlgebra(transformation_matrix).diagonalize
        orbital_energies, orbital_coefficients = diagonalize(core_hamiltonian)

        print('\nORBITAL ENERGY EIGENVALUES')
        print(orbital_energies)
        print('\nORBITAL COEFFICIENTS')
        print(orbital_coefficients)

        print('\n*****************************************************************************************************')
        print('\nBEGIN TWO ELECTRON REPULSION CALCULATION')
        start_repulsion = time.clock()
        # repulsion_dictionary = TwoElectronRepulsionElementCook(basis_set_array).store_parallel(4)
        # repulsion_dictionary = TwoElectronRepulsionElementHGP(basis_set_array).store_parallel(4)
        repulsion = TwoElectronRepulsionElementOS(self.basis_set_array).store_parallel(4)
        print('TIME TAKEN: ' + str(time.clock() - start_repulsion) + 's\n')

        scf = SelfConsistentField(core_hamiltonian, repulsion, diagonalize, create_matrix, self.electrons, self.multiplicity)
        return scf, orbital_coefficients, repulsion

    def restricted(self):
        scf, coefficients_initial_guess, repulsion = self.begin()

        print('\nBEGIN SCF PROCEDURE')
        start = time.clock()
        electron_energy, orbital_energies, orbital_coefficients = scf.restricted(coefficients_initial_guess)
        print('TIME TAKEN: ' + str(time.clock() - start) + 's\n')

        print('\nORBITAL ENERGY EIGENVALUES')
        print(orbital_energies)
        print('\nORBITAL COEFFICIENTS')
        print(orbital_coefficients, end='\n\n')
        return electron_energy, orbital_energies, orbital_coefficients, repulsion

    def unrestricted(self):
        scf, coefficients_initial_guess, repulsion = self.begin()

        print('\nBEGIN SCF PROCEDURE')
        start = time.clock()
        electron_energy, energies_alpha, energies_beta, coefficients_alpha, coefficients_beta = scf.unrestricted(coefficients_initial_guess)
        print('TIME TAKEN: ' + str(time.clock() - start) + 's\n')

        print('\nALPHA ORBITAL ENERGY EIGENVALUES')
        print(energies_alpha)
        print('\nBETA ORBITAL ENERGY EIGENVALUES')
        print(energies_beta)
        print('\nALPHA ORBITAL COEFFICIENTS')
        print(coefficients_alpha, end='\n\n')
        print('\nBETA ORBITAL COEFFICIENTS')
        print(coefficients_beta, end='\n\n')
        return electron_energy, energies_alpha, energies_beta, coefficients_alpha, coefficients_beta, repulsion
