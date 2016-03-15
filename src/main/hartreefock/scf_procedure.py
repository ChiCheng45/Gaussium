from src.main.matrixelements import DensityMatrixRestricted
from src.main.matrixelements import DensityMatrixUnrestricted
from src.main.matrixelements import GMatrixRestricted
from src.main.matrixelements import GMatrixUnrestricted
from src.main.hartreefock import TotalEnergy
from src.main.diis.diis import DIIS
from math import floor, ceil


class SelfConsistentField:

    def __init__(self, core_hamiltonian, linear_algebra, density_matrix_factory, g_matrix_factory, diis):
        self.core_hamiltonian = core_hamiltonian
        self.linear_algebra = linear_algebra
        self.density_matrix_factory = density_matrix_factory
        self.g_matrix_factory = g_matrix_factory
        self.diis = diis
        self.calculate = TotalEnergy(core_hamiltonian)
        self.total_energy = 0
        self.previous_total_energy = 0
        self.delta_energy = 1


class RestrictedSCF(SelfConsistentField):

    def __init__(self, core_hamiltonian, linear_algebra, repulsion, electrons, overlap):
        super().__init__(core_hamiltonian, linear_algebra, DensityMatrixRestricted(electrons), GMatrixRestricted(repulsion), DIIS(overlap))

    def begin(self, orbital_coefficients):
        orbital_energies = []

        while abs(self.delta_energy) > 1e-9:
            density_matrix = self.density_matrix_factory.create(orbital_coefficients)
            g_matrix = self.g_matrix_factory.create(density_matrix)
            fock_matrix = self.core_hamiltonian + g_matrix
            self.total_energy = self.calculate.restricted(density_matrix, fock_matrix)
            self.delta_energy = self.previous_total_energy - self.total_energy
            self.previous_total_energy = self.total_energy
            print('SCF ENERGY: ' + str(self.total_energy) + ' a.u.')

            if abs(self.delta_energy) > 1e-9:
                fock_matrix = self.diis.fock_matrix(fock_matrix, density_matrix)
                orbital_energies, orbital_coefficients = self.linear_algebra.diagonalize(fock_matrix)

        return self.total_energy, orbital_energies, orbital_coefficients


class UnrestrictedSCF(SelfConsistentField):

    def __init__(self, core_hamiltonian, linear_algebra, repulsion, electrons, multiplicity):
        super().__init__(core_hamiltonian, linear_algebra, DensityMatrixUnrestricted(), GMatrixUnrestricted(repulsion))
        self.electrons = electrons
        self.multiplicity = multiplicity

    def begin(self, orbital_coefficients):
        difference = floor(self.multiplicity / 2)
        electrons_alph = ceil(self.electrons / 2) + difference
        electrons_beta = floor(self.electrons / 2) - difference
        coefficients_alph = orbital_coefficients
        coefficients_beta = orbital_coefficients
        energies_alph = []
        energies_beta = []

        while abs(self.delta_energy) > 1e-9:

            if self.total_energy == 0:
                density_matrix_alph = self.density_matrix_factory.create((electrons_alph + 1), coefficients_alph)
                density_matrix_beta = self.density_matrix_factory.create((electrons_beta - 1), coefficients_beta)
            else:
                density_matrix_alph = self.density_matrix_factory.create(electrons_alph, coefficients_alph)
                density_matrix_beta = self.density_matrix_factory.create(electrons_beta, coefficients_beta)

            self.g_matrix_factory.set_density_matrices(density_matrix_alph, density_matrix_beta)
            g_matrix_alph = self.g_matrix_factory.create_alph()
            g_matrix_beta = self.g_matrix_factory.create_beta()
            fock_matrix_alph = self.core_hamiltonian + g_matrix_alph
            fock_matrix_beta = self.core_hamiltonian + g_matrix_beta
            energies_alph, coefficients_alph = self.linear_algebra.diagonalize(fock_matrix_alph)
            energies_beta, coefficients_beta = self.linear_algebra.diagonalize(fock_matrix_beta)

            self.total_energy = self.calculate.unrestricted(density_matrix_alph, density_matrix_beta, fock_matrix_alph, fock_matrix_beta)
            self.delta_energy = self.previous_total_energy - self.total_energy
            self.previous_total_energy = self.total_energy
            print('SCF ENERGY: ' + str(self.total_energy) + ' a.u.')

        return self.total_energy, energies_alph, energies_beta, coefficients_alph, coefficients_beta
