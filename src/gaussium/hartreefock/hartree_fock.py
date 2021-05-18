import time
import numpy as np
from gaussium.hartreefock import BlockedFockMatrixUnrestricted
from gaussium.hartreefock import BlockedLinearAlgebra
from gaussium.hartreefock import BlockedUnrestrictedSCF
from gaussium.hartreefock import FockMatrixRestricted
from gaussium.hartreefock import FockMatrixUnrestricted
from gaussium.hartreefock import LinearAlgebra
from gaussium.hartreefock import PopleNesbetBerthier
from gaussium.hartreefock import RestrictedSCF
from gaussium.matrixelements import KineticEnergyMatrix
from gaussium.matrixelements import NuclearAttractionMatrix
from gaussium.matrixelements import OrbitalOverlapMatrix
from gaussium.matrixelements import TwoElectronRepulsionMatrixOS
from gaussium.matrixelements import blocked_spin_basis_set


class HartreeFock:

    def __init__(self, nuclei_array, basis_set_array, electrons, symmetry, processes):
        self.scf_method = None
        self.nuclei_array = nuclei_array
        self.basis_set_array = basis_set_array
        self.electrons = electrons
        self.symmetry = symmetry
        self.orbital_overlap = OrbitalOverlapMatrix(basis_set_array).create()
        self.kinetic_energy = KineticEnergyMatrix(basis_set_array).create()
        self.nuclear_attraction = NuclearAttractionMatrix(basis_set_array, nuclei_array).create()
        self.core_hamiltonian = self.kinetic_energy + self.nuclear_attraction
        self.linear_algebra = LinearAlgebra(self.orbital_overlap)
        print('\n*************************************************************************************************')
        print('\nMATRICES\n')
        print('\nORBITAL OVERLAP MATRIX\n{}'.format(self.orbital_overlap))
        print('\nKINETIC ENERGY MATRIX\n{}'.format(self.kinetic_energy))
        print('\nNUCLEAR POTENTIAL ENERGY MATRIX\n{}'.format(self.nuclear_attraction))
        print('\nCORE HAMILTONIAN MATRIX\n{}'.format(self.core_hamiltonian))
        print('\nBEGIN TWO ELECTRON REPULSION CALCULATION')
        start_repulsion = time.perf_counter()
        self.repulsion = TwoElectronRepulsionMatrixOS(
            self.basis_set_array, self.symmetry, processes
        ).create_repulsion_matrix()
        print('TIME TAKEN: ' + str(time.perf_counter() - start_repulsion) + 's\n')
        print('\n*************************************************************************************************')

    def initial_guess(self):
        initial_orbital_energies, initial_orbital_coefficients = self.linear_algebra.diagonalize(self.core_hamiltonian)
        return initial_orbital_energies, initial_orbital_coefficients

    def begin_scf(self):
        pass

    def energies(self):
        energy = self.begin_scf()[0]
        return energy, 0.0


class Restricted(HartreeFock):

    def __init__(self, nuclei_array, basis_set_array, electrons, symmetry, processes):
        super().__init__(nuclei_array, basis_set_array, electrons, symmetry, processes)

    def begin_scf(self):
        initial_energies, initial_coefficients = self.initial_guess()
        print('COEFFICIENTS INITIAL GUESS\n{}'.format(initial_coefficients))
        print('\n\nBEGIN SCF PROCEDURE')
        start = time.perf_counter()
        electron_energy, orbital_energies, orbital_coefficients \
            = self.scf_method.begin_iterations(initial_energies, initial_coefficients)
        print('TIME TAKEN: ' + str(time.perf_counter() - start) + 's\n')
        print('\nORBITAL ENERGY EIGENVALUES\n{}'.format(orbital_energies))
        print('\nORBITAL COEFFICIENTS\n{}'.format(orbital_coefficients), end='\n\n\n')

        return electron_energy, orbital_energies, orbital_coefficients


class RestrictedHF(Restricted):

    def __init__(self, nuclei_array, basis_set_array, electrons, symmetry, processes):
        super().__init__(nuclei_array, basis_set_array, electrons, symmetry, processes)
        self.scf_method = RestrictedSCF(
            self.linear_algebra, self.electrons, self.orbital_overlap,
            FockMatrixRestricted(self.core_hamiltonian, self.repulsion)
        )
        print('\nBEGIN RESTRICTED HARTREE FOCK\n')


class Unrestricted(HartreeFock):

    def __init__(self, nuclei_array, basis_set_array, electrons, symmetry, processes):
        super().__init__(nuclei_array, basis_set_array, electrons, symmetry, processes)

    def begin_scf(self):
        initial_energies, initial_coefficients = self.initial_guess()
        print('COEFFICIENTS INITIAL GUESS\n{}'.format(initial_coefficients))
        print('\n\nBEGIN SCF PROCEDURE')
        start = time.perf_counter()
        electron_energy, energies_alpha, energies_beta, coefficients_alpha, coefficients_beta \
            = self.scf_method.begin_iterations(initial_energies, initial_coefficients)
        print('TIME TAKEN: ' + str(time.perf_counter() - start) + 's\n')
        print('\nALPHA ORBITAL ENERGY EIGENVALUES\n{}'.format(energies_alpha))
        print('\nBETA ORBITAL ENERGY EIGENVALUES\n{}'.format(energies_beta))
        print('\nALPHA ORBITAL COEFFICIENTS\n{}'.format(coefficients_alpha), end='\n')
        print('\nBETA ORBITAL COEFFICIENTS\n{}'.format(coefficients_beta), end='\n\n\n')

        return electron_energy, energies_alpha, energies_beta, coefficients_alpha, coefficients_beta


class UnrestrictedHF(Unrestricted):

    def __init__(self, nuclei_array, basis_set_array, electrons, multiplicity, symmetry, processes):
        super().__init__(nuclei_array, basis_set_array, electrons, symmetry, processes)
        self.scf_method = PopleNesbetBerthier(
            self.linear_algebra, self.electrons, multiplicity,
            FockMatrixUnrestricted(self.core_hamiltonian, self.repulsion)
        )
        print('\nBEGIN UNRESTRICTED HARTREE FOCK\n')


class BlockedHartreeFock(Restricted):

    def __init__(self, nuclei_array, basis_set_array, electrons, multiplicity, symmetry, processes):
        super().__init__(nuclei_array, basis_set_array, electrons, symmetry, processes)
        self.zeros = np.zeros((self.orbital_overlap.shape[0], self.orbital_overlap.shape[0]))

        self.orbital_overlap = np.block([
                [self.orbital_overlap, self.zeros],
                [self.zeros, self.orbital_overlap]
        ])

        self.core_hamiltonian = np.block([
                [self.core_hamiltonian, self.zeros],
                [self.zeros, self.core_hamiltonian]
        ])

        self.repulsion = blocked_spin_basis_set(self.repulsion)
        self.linear_algebra = BlockedLinearAlgebra(self.orbital_overlap)
        self.scf_method = BlockedUnrestrictedSCF(self.linear_algebra, self.electrons, multiplicity,
        self.orbital_overlap, BlockedFockMatrixUnrestricted(self.core_hamiltonian, self.repulsion))
        print('\nBEGIN BLOCKED UNRESTRICTED HARTREE FOCK\n')
