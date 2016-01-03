from src.main.matrixelements import DensityMatrixRestricted, DensityMatrixUnrestricted
from src.main.matrixelements import GMatrixRestricted, GMatrixUnrestricted
from src.main.matrixelements import KineticEnergyMatrix
from src.main.matrixelements import NuclearAttractionMatrix
from src.main.matrixelements import OrbitalOverlapMatrix
from src.main.matrixelements import TwoElectronRepulsionElementCook
from src.main.matrixelements import TwoElectronRepulsionElementHGP
from src.main.matrixelements import TwoElectronRepulsionElementOS
from src.main.hartreefock import SelfConsistentField, LinearAlgebra
import numpy as np
import time


class HartreeFock:

    def __init__(self, nuclei_array, electrons, multiplicity, basis_set_array):
        self.nuclei_array = nuclei_array
        self.electrons = electrons
        self.multiplicity = multiplicity
        self.basis_set_array = basis_set_array
        self.core_hamiltonian = np.matrix([])
        self.repulsion = {}

    def begin(self):
        print('\n*****************************************************************************************************')
        orbital_overlap = OrbitalOverlapMatrix().create(self.basis_set_array)
        kinetic_energy = KineticEnergyMatrix().create(self.basis_set_array)
        nuclear_potential = NuclearAttractionMatrix().create(self.basis_set_array, self.nuclei_array)
        self.core_hamiltonian = kinetic_energy + nuclear_potential
        transformation_matrix = LinearAlgebra.transformation_matrix(orbital_overlap)

        print('\nMATRICES\n')
        print('\nORBITAL OVERLAP MATRIX')
        print(orbital_overlap)
        print('\nKINETIC ENERGY MATRIX')
        print(kinetic_energy)
        print('\nNUCLEAR POTENTIAL ENERGY MATRIX')
        print(nuclear_potential)
        print('\nCORE HAMILTONIAN MATRIX')
        print(self.core_hamiltonian)
        print('\nTRANSFORMATION MATRIX')
        print(transformation_matrix)

        print('\n*****************************************************************************************************')
        print('\nINITIAL GUESS\n')
        linear_algebra = LinearAlgebra(transformation_matrix)
        initial_orbital_energies, initial_orbital_coefficients = linear_algebra.diagonalize(self.core_hamiltonian)

        print('\nORBITAL ENERGY EIGENVALUES')
        print(initial_orbital_energies)
        print('\nORBITAL COEFFICIENTS')
        print(initial_orbital_coefficients)

        print('\n*****************************************************************************************************')
        print('\nBEGIN TWO ELECTRON REPULSION CALCULATION')
        start_repulsion = time.clock()
        # repulsion_dictionary = TwoElectronRepulsionElementCook(basis_set_array).store_parallel(4)
        # repulsion_dictionary = TwoElectronRepulsionElementHGP(basis_set_array).store_parallel(4)
        self.repulsion = TwoElectronRepulsionElementOS(self.basis_set_array).store_parallel(4)
        print('TIME TAKEN: ' + str(time.clock() - start_repulsion) + 's\n')

        return initial_orbital_coefficients, linear_algebra

    def restricted(self):
        initial_orbital_coefficients, linear_algebra = self.begin()
        scf = SelfConsistentField(self.core_hamiltonian, DensityMatrixRestricted(), GMatrixRestricted(self.repulsion), linear_algebra, self.electrons, self.multiplicity)

        print('\nBEGIN SCF PROCEDURE')
        start = time.clock()
        electron_energy, orbital_energies, orbital_coefficients = scf.restricted(initial_orbital_coefficients)
        print('TIME TAKEN: ' + str(time.clock() - start) + 's\n')

        print('\nORBITAL ENERGY EIGENVALUES')
        print(orbital_energies)
        print('\nORBITAL COEFFICIENTS')
        print(orbital_coefficients, end='\n\n')
        return electron_energy, orbital_energies, orbital_coefficients, self.repulsion

    def unrestricted(self):
        initial_orbital_coefficients, linear_algebra = self.begin()
        scf = SelfConsistentField(self.core_hamiltonian, DensityMatrixUnrestricted(), GMatrixUnrestricted(self.repulsion), linear_algebra, self.electrons, self.multiplicity)

        print('\nBEGIN SCF PROCEDURE')
        start = time.clock()
        electron_energy, energies_alpha, energies_beta, coefficients_alpha, coefficients_beta = scf.unrestricted(initial_orbital_coefficients)
        print('TIME TAKEN: ' + str(time.clock() - start) + 's\n')

        print('\nALPHA ORBITAL ENERGY EIGENVALUES')
        print(energies_alpha)
        print('\nBETA ORBITAL ENERGY EIGENVALUES')
        print(energies_beta)
        print('\nALPHA ORBITAL COEFFICIENTS')
        print(coefficients_alpha, end='\n')
        print('\nBETA ORBITAL COEFFICIENTS')
        print(coefficients_beta, end='\n\n')
        return electron_energy, energies_alpha, energies_beta, coefficients_alpha, coefficients_beta, self.repulsion
