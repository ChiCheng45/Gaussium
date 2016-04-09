from src.main.matrixelements import DensityMatrixRestricted
from src.main.matrixelements import DensityMatrixUnrestricted
from src.main.matrixelements import FockMatrixRestricted
from src.main.matrixelements import FockMatrixUnrestricted
from src.main.matrixelements import FockMatrixConstrained
from src.main.hartreefock import TotalEnergy
from src.main.diismethod import DIIS
from math import floor, ceil


class SelfConsistentField:

    def __init__(self, core_hamiltonian, linear_algebra, electrons, density_matrix_factory, fock_matrix_factory):
        self.linear_algebra = linear_algebra
        self.electrons = electrons
        self.density_matrix_factory = density_matrix_factory
        self.fock_matrix_factory = fock_matrix_factory
        self.calculate = TotalEnergy(core_hamiltonian)
        self.total_energy = 0
        self.previous_total_energy = 0
        self.delta_energy = 10


class RestrictedSCF(SelfConsistentField):

    def __init__(self, core_hamiltonian, linear_algebra, repulsion, electrons, overlap):
        super().__init__(core_hamiltonian, linear_algebra, electrons, DensityMatrixRestricted(electrons), FockMatrixRestricted(core_hamiltonian, repulsion))
        self.diis = DIIS(overlap, linear_algebra)

    def begin(self, orbital_coefficients):
        orbital_energies = []

        while abs(self.delta_energy) > 1e-12:
            density_matrix = self.density_matrix_factory.create(orbital_coefficients)
            fock_matrix = self.fock_matrix_factory.create(density_matrix)
            self.total_energy = self.calculate.restricted(density_matrix, fock_matrix)
            self.delta_energy = self.previous_total_energy - self.total_energy
            self.previous_total_energy = self.total_energy
            print('SCF ENERGY: ' + str(self.total_energy) + ' a.u.')

            if abs(self.delta_energy) > 1e-12:
                fock_matrix = self.diis.fock_matrix(fock_matrix, density_matrix)
                orbital_energies, orbital_coefficients = self.linear_algebra.diagonalize(fock_matrix)

        return self.total_energy, orbital_energies, orbital_coefficients


class PopleNesbetBerthier(SelfConsistentField):

    def __init__(self, core_hamiltonian, linear_algebra, electrons, multiplicity, density_matrix_factory, fock_matrix_factory, overlap):
        super().__init__(core_hamiltonian, linear_algebra, electrons, density_matrix_factory, fock_matrix_factory)
        self.electrons_alph = ceil(self.electrons / 2) + floor(multiplicity / 2)
        self.electrons_beta = floor(self.electrons / 2) - floor(multiplicity / 2)
        self.diis_alph = DIIS(overlap, linear_algebra)
        self.diis_beta = DIIS(overlap, linear_algebra)

    def begin(self, orbital_coefficients):
        coefficients_alph = orbital_coefficients
        coefficients_beta = orbital_coefficients
        energies_alph = []
        energies_beta = []

        while abs(self.delta_energy) > 1e-12:

            if self.total_energy == 0:
                density_matrix_alph = self.density_matrix_factory.create((self.electrons_alph + 1), coefficients_alph)
                density_matrix_beta = self.density_matrix_factory.create((self.electrons_beta - 1), coefficients_beta)
            else:
                density_matrix_alph = self.density_matrix_factory.create(self.electrons_alph, coefficients_alph)
                density_matrix_beta = self.density_matrix_factory.create(self.electrons_beta, coefficients_beta)

            fock_matrix_alph, fock_matrix_beta = self.fock_matrix_factory.create(density_matrix_alph, density_matrix_beta)

            self.total_energy = self.calculate.unrestricted(density_matrix_alph, density_matrix_beta, fock_matrix_alph, fock_matrix_beta)
            self.delta_energy = self.previous_total_energy - self.total_energy
            self.previous_total_energy = self.total_energy
            print('SCF ENERGY: ' + str(self.total_energy) + ' a.u.')

            if abs(self.delta_energy) > 1e-12:
                fock_matrix_alph = self.diis_alph.fock_matrix(fock_matrix_alph, density_matrix_alph)
                fock_matrix_beta = self.diis_beta.fock_matrix(fock_matrix_beta, density_matrix_beta)
                energies_alph, coefficients_alph = self.linear_algebra.diagonalize(fock_matrix_alph)
                energies_beta, coefficients_beta = self.linear_algebra.diagonalize(fock_matrix_beta)

        return self.total_energy, energies_alph, energies_beta, coefficients_alph, coefficients_beta


class DifferentOrbitalsDifferentSpins(PopleNesbetBerthier):

    def __init__(self, core_hamiltonian, linear_algebra, repulsion, electrons, multiplicity, overlap):
        super().__init__(core_hamiltonian, linear_algebra, electrons, multiplicity, DensityMatrixUnrestricted(),
                         FockMatrixUnrestricted(core_hamiltonian, repulsion), overlap)


class ConstrainedUnrestrictedSCF(PopleNesbetBerthier):

    def __init__(self, core_hamiltonian, linear_algebra, repulsion, electrons, multiplicity, overlap):
        super().__init__(core_hamiltonian, linear_algebra, electrons, multiplicity, DensityMatrixUnrestricted(),
                         FockMatrixConstrained(core_hamiltonian, repulsion, electrons, multiplicity, linear_algebra),
                         overlap)
