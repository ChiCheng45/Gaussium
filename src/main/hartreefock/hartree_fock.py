from src.main.matrixelements import KineticEnergyMatrix
from src.main.matrixelements import NuclearAttractionMatrix
from src.main.matrixelements import OrbitalOverlapMatrix
from src.main.matrixelements import TwoElectronRepulsionMatrixOS
from src.main.hartreefock import LinearAlgebra
from src.main.hartreefock import BlockedLinearAlgebra
from src.main.hartreefock import RestrictedSCF
from src.main.hartreefock import DifferentOrbitalsDifferentSpins
from src.main.hartreefock import ConstrainedUnrestrictedSCF
from src.main.hartreefock import BlockedUnrestrictedSCF
from src.main.matrixelements import IntegralTransformations
import numpy as np
import time


class HartreeFock:

    def __init__(self, nuclei_array, basis_set_array, electrons, scf_method):
        self.nuclei_array = nuclei_array
        self.basis_set_array = basis_set_array
        self.electrons = electrons
        self.scf_method = scf_method
        self.orbital_overlap_matrix = OrbitalOverlapMatrix(basis_set_array)
        self.kinetic_energy_matrix = KineticEnergyMatrix(basis_set_array)
        self.nuclear_attraction_matrix = NuclearAttractionMatrix(basis_set_array, nuclei_array)
        self.repulsion_elements = TwoElectronRepulsionMatrixOS(basis_set_array)
        self.linear_algebra = LinearAlgebra
        self.core_hamiltonian = np.matrix([])
        self.orbital_overlap = np.matrix([])
        self.repulsion = np.matrix([])

    def start(self):
        print('\n*****************************************************************************************************')
        kinetic_energy = self.kinetic_energy_matrix.create()
        nuclear_potential = self.nuclear_attraction_matrix.create()
        self.orbital_overlap = self.orbital_overlap_matrix.create()
        self.linear_algebra = self.linear_algebra(self.orbital_overlap)
        self.core_hamiltonian = kinetic_energy + nuclear_potential

        print('\nMATRICES\n')
        print('\nORBITAL OVERLAP MATRIX')
        print(self.orbital_overlap)
        print('\nKINETIC ENERGY MATRIX')
        print(kinetic_energy)
        print('\nNUCLEAR POTENTIAL ENERGY MATRIX')
        print(nuclear_potential)
        print('\nCORE HAMILTONIAN MATRIX')
        print(self.core_hamiltonian)

        print('\n*****************************************************************************************************')
        print('\nINITIAL GUESS\n')
        initial_orbital_energies, initial_orbital_coefficients = self.linear_algebra.diagonalize(self.core_hamiltonian)

        print('\nORBITAL ENERGY EIGENVALUES')
        print(initial_orbital_energies)
        print('\nORBITAL COEFFICIENTS')
        print(initial_orbital_coefficients)

        print('\n*****************************************************************************************************')
        print('\nBEGIN TWO ELECTRON REPULSION CALCULATION')
        start_repulsion = time.clock()
        self.repulsion = self.repulsion_elements.store_parallel(4)
        print('TIME TAKEN: ' + str(time.clock() - start_repulsion) + 's\n')

        return initial_orbital_coefficients


class RestrictedHF(HartreeFock):

    def __init__(self, nuclei_array, basis_set_array, electrons):
        super().__init__(nuclei_array, basis_set_array, electrons, RestrictedSCF)

    def begin(self):
        initial_coefficients = self.start()
        self.scf_method = self.scf_method(self.core_hamiltonian, self.linear_algebra, self.repulsion, self.electrons,
        self.orbital_overlap)

        start = time.clock()
        print('\nBEGIN SCF PROCEDURE')
        electron_energy, orbital_energies, orbital_coefficients = self.scf_method.begin(initial_coefficients)
        print('TIME TAKEN: ' + str(time.clock() - start) + 's\n')

        print('\nORBITAL ENERGY EIGENVALUES')
        print(orbital_energies)
        print('\nORBITAL COEFFICIENTS')
        print(orbital_coefficients, end='\n\n\n')
        return electron_energy, orbital_energies, orbital_coefficients, self.repulsion


class UnrestrictedHF(HartreeFock):

    def __init__(self, nuclei_array, basis_set_array, electrons, multiplicity, scf_method):
        super().__init__(nuclei_array, basis_set_array, electrons, scf_method)
        self.multiplicity = multiplicity

    def begin(self):
        initial_coefficients = self.start()
        self.scf_method = self.scf_method(self.core_hamiltonian, self.linear_algebra, self.repulsion, self.electrons,
        self.multiplicity, self.orbital_overlap)

        start = time.clock()
        print('\nBEGIN SCF PROCEDURE')
        electron_energy, energies_alpha, energies_beta, coefficients_alpha, coefficients_beta \
        = self.scf_method.begin(initial_coefficients)
        print('TIME TAKEN: ' + str(time.clock() - start) + 's\n')

        print('\nALPHA ORBITAL ENERGY EIGENVALUES')
        print(energies_alpha)
        print('\nBETA ORBITAL ENERGY EIGENVALUES')
        print(energies_beta)
        print('\nALPHA ORBITAL COEFFICIENTS')
        print(coefficients_alpha, end='\n')
        print('\nBETA ORBITAL COEFFICIENTS')
        print(coefficients_beta, end='\n\n\n')
        return electron_energy, energies_alpha, energies_beta, coefficients_alpha, coefficients_beta, self.repulsion


class DODSUnrestricted(UnrestrictedHF):

    def __init__(self, nuclei_array, basis_set_array, electrons, multiplicity):
        super().__init__(nuclei_array, basis_set_array, electrons, multiplicity, DifferentOrbitalsDifferentSpins)


class ConstrainedUnrestricted(UnrestrictedHF):

    def __init__(self, nuclei_array, basis_set_array, electrons, multiplicity):
        super().__init__(nuclei_array, basis_set_array, electrons, multiplicity, ConstrainedUnrestrictedSCF)


class BlockedHartreeFock(HartreeFock):

    def __init__(self, nuclei_array, basis_set_array, electrons, multiplicity):
        super().__init__(nuclei_array, basis_set_array, electrons, BlockedUnrestrictedSCF)
        self.multiplicity = multiplicity
        self.block_linear_algebra = BlockedLinearAlgebra

    def begin(self):
        initial_coefficients = self.start()
        zeros = np.zeros((self.orbital_overlap.shape[0], self.orbital_overlap.shape[0]))

        self.orbital_overlap = np.bmat([
                [self.orbital_overlap, zeros],
                [zeros, self.orbital_overlap]])

        self.core_hamiltonian = np.bmat([
                [self.core_hamiltonian, zeros],
                [zeros, self.core_hamiltonian]])

        initial_coefficients = np.bmat([
                [initial_coefficients, zeros],
                [zeros, initial_coefficients]])

        self.repulsion = IntegralTransformations.create_spin_basis_integrals(self.repulsion)
        self.block_linear_algebra = self.block_linear_algebra(self.orbital_overlap)

        self.scf_method = self.scf_method(self.core_hamiltonian, self.block_linear_algebra, self.repulsion,
        self.electrons, self.multiplicity, self.orbital_overlap)

        start = time.clock()
        print('\nBEGIN SCF PROCEDURE')
        electron_energy, orbital_energies, orbital_coefficients = self.scf_method.begin(initial_coefficients)
        print('TIME TAKEN: ' + str(time.clock() - start) + 's\n')

        print('\nORBITAL ENERGY EIGENVALUES')
        print(orbital_energies)
        print('\nORBITAL COEFFICIENTS')
        print(orbital_coefficients, end='\n\n\n')
        return electron_energy, orbital_energies, orbital_coefficients, self.repulsion
