from src.main.matrixelements import density_matrix_unrestricted
from src.main.matrixelements import density_matrix_restricted
from src.main.matrixelements import blocked_density_matrix
from src.main.matrixelements import FockMatrixRestricted
from src.main.matrixelements import FockMatrixUnrestricted
from src.main.matrixelements import FockMatrixConstrained
from src.main.matrixelements import BlockedFockMatrixUnrestricted
from src.main.hartreefock import TotalEnergy
from src.main.diismethod import DIIS


class SelfConsistentField:

    def __init__(self, core_hamiltonian, linear_algebra, electrons, fock_matrix_factory, threshold=1e-12):
        self.linear_algebra = linear_algebra
        self.electrons = electrons
        self.fock_matrix_factory = fock_matrix_factory
        self.calculate = TotalEnergy(core_hamiltonian)
        self.threshold = threshold


class RestrictedSCF(SelfConsistentField):

    def __init__(self, core_hamiltonian, linear_algebra, repulsion, electrons, overlap):
        super().__init__(core_hamiltonian, linear_algebra, electrons, FockMatrixRestricted(core_hamiltonian, repulsion))
        self.diis = DIIS(overlap, linear_algebra)

    def begin_iterations(self, orbital_coefficients):
        orbital_energies = []
        total_energy = previous_total_energy = 0

        while True:

            density_matrix = density_matrix_restricted(orbital_coefficients, self.electrons)
            fock_matrix = self.fock_matrix_factory.create(density_matrix)
            total_energy = self.calculate.restricted(density_matrix, fock_matrix)
            delta_energy = previous_total_energy - total_energy
            previous_total_energy = total_energy
            print('SCF ENERGY: ' + str(total_energy) + ' a.u.')

            if abs(delta_energy) < self.threshold:
                break

            fock_matrix = self.diis.fock_matrix(fock_matrix, density_matrix)
            orbital_energies, orbital_coefficients = self.linear_algebra.diagonalize(fock_matrix)

        return total_energy, orbital_energies, orbital_coefficients


class PopleNesbetBerthier(SelfConsistentField):

    def __init__(self, core_hamiltonian, linear_algebra, electrons, multiplicity, fock_matrix_factory):
        super().__init__(core_hamiltonian, linear_algebra, electrons, fock_matrix_factory)
        self.electrons_alph = (electrons + multiplicity - 1) // 2
        self.electrons_beta = (electrons - multiplicity + 1) // 2

    def begin_iterations(self, orbital_coefficients):
        coefficients_alph = orbital_coefficients
        coefficients_beta = orbital_coefficients
        energies_alph = []
        energies_beta = []
        total_energy = previous_total_energy = 0

        while True:

            if total_energy == 0:
                density_matrix_alph = density_matrix_unrestricted(coefficients_alph, coefficients_alph.shape[0])
                density_matrix_beta = density_matrix_unrestricted(coefficients_beta, 0)
            else:
                density_matrix_alph = density_matrix_unrestricted(coefficients_alph, self.electrons_alph)
                density_matrix_beta = density_matrix_unrestricted(coefficients_beta, self.electrons_beta)

            fock_matrix_alph, fock_matrix_beta = self.fock_matrix_factory.create(density_matrix_alph,
            density_matrix_beta)

            total_energy = self.calculate.unrestricted(density_matrix_alph, density_matrix_beta, fock_matrix_alph,
            fock_matrix_beta)
            delta_energy = previous_total_energy - total_energy
            previous_total_energy = total_energy
            print('SCF ENERGY: ' + str(total_energy) + ' a.u.')

            if abs(delta_energy) < self.threshold:
                break

            energies_alph, coefficients_alph = self.linear_algebra.diagonalize(fock_matrix_alph)
            energies_beta, coefficients_beta = self.linear_algebra.diagonalize(fock_matrix_beta)

        return total_energy, energies_alph, energies_beta, coefficients_alph, coefficients_beta


class DifferentOrbitalsDifferentSpins(PopleNesbetBerthier):

    def __init__(self, core_hamiltonian, linear_algebra, repulsion, electrons, multiplicity):
        super().__init__(core_hamiltonian, linear_algebra, electrons, multiplicity,
        FockMatrixUnrestricted(core_hamiltonian, repulsion))


class ConstrainedUnrestrictedSCF(PopleNesbetBerthier):

    def __init__(self, core_hamiltonian, linear_algebra, repulsion, electrons, multiplicity):
        super().__init__(core_hamiltonian, linear_algebra, electrons, multiplicity,
        FockMatrixConstrained(core_hamiltonian, repulsion, electrons, multiplicity, linear_algebra))


class BlockedUnrestrictedSCF(SelfConsistentField):

    def __init__(self, core_hamiltonian, linear_algebra, repulsion, electrons, multiplicity, overlap):
        super().__init__(core_hamiltonian, linear_algebra, electrons,
        BlockedFockMatrixUnrestricted(core_hamiltonian, repulsion))
        self.electrons_alph = (electrons + multiplicity - 1) // 2
        self.electrons_beta = (electrons - multiplicity + 1) // 2
        self.diis = DIIS(overlap, linear_algebra)

    def begin_iterations(self, orbital_coefficients):
        orbital_energies = []
        total_energy = previous_total_energy = 0

        while True:

            if total_energy == 0:
                density_matrix = blocked_density_matrix(orbital_coefficients, orbital_coefficients.shape[0] // 2, 0)
            else:
                density_matrix = blocked_density_matrix(orbital_coefficients, self.electrons_alph, self.electrons_beta)

            fock_matrix = self.fock_matrix_factory.create(density_matrix)
            total_energy = self.calculate.restricted(density_matrix, fock_matrix)
            delta_energy = previous_total_energy - total_energy
            previous_total_energy = total_energy
            print('SCF ENERGY: ' + str(total_energy) + ' a.u.')

            if abs(delta_energy) < self.threshold:
                break

            fock_matrix = self.diis.fock_matrix(fock_matrix, density_matrix)
            orbital_energies, orbital_coefficients = self.linear_algebra.diagonalize(fock_matrix)

        return total_energy, orbital_energies, orbital_coefficients
