from src.diismethod import DIIS
from src.hartreefock import TotalEnergy
from src.matrixelements import blocked_density_matrix
from src.matrixelements import density_matrix_restricted
from src.matrixelements import density_matrix_unrestricted


class SelfConsistentField:

    def __init__(self, linear_algebra, electrons, hamiltonian_matrix_factory, threshold=1e-12):
        self.linear_algebra = linear_algebra
        self.electrons = electrons
        self.hamiltonian_matrix_factory = hamiltonian_matrix_factory
        self.threshold = threshold


class RestrictedSCF(SelfConsistentField):

    def __init__(self, linear_algebra, electrons, overlap, hamiltonian_matrix_factory):
        super().__init__(linear_algebra, electrons, hamiltonian_matrix_factory)
        self.calculate = TotalEnergy(hamiltonian_matrix_factory.core_hamiltonian)
        self.diis = DIIS(overlap, linear_algebra)

    def begin_iterations(self, orbital_energies, orbital_coefficients):
        previous_total_energy = 0

        while True:

            density_matrix = density_matrix_restricted(orbital_coefficients, self.electrons)
            fock_matrix = self.hamiltonian_matrix_factory.create(density_matrix)
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

    def __init__(self, linear_algebra, electrons, multiplicity, hamiltonian_matrix_factory):
        super().__init__(linear_algebra, electrons, hamiltonian_matrix_factory)
        self.calculate = TotalEnergy(hamiltonian_matrix_factory.core_hamiltonian)
        self.electrons_alph = (electrons + multiplicity - 1) // 2
        self.electrons_beta = (electrons - multiplicity + 1) // 2

    def begin_iterations(self, orbital_energies, orbital_coefficients):
        coefficients_alph = orbital_coefficients
        coefficients_beta = orbital_coefficients
        energies_alph = orbital_energies
        energies_beta = orbital_energies
        total_energy = previous_total_energy = 0

        while True:

            if total_energy == 0:
                density_matrix_alph = density_matrix_unrestricted(coefficients_alph, coefficients_alph.shape[0])
                density_matrix_beta = density_matrix_unrestricted(coefficients_beta, 0)
            else:
                density_matrix_alph = density_matrix_unrestricted(coefficients_alph, self.electrons_alph)
                density_matrix_beta = density_matrix_unrestricted(coefficients_beta, self.electrons_beta)

            fock_matrix_alph, fock_matrix_beta = self.hamiltonian_matrix_factory.create(
                density_matrix_alph, density_matrix_beta
            )

            total_energy = self.calculate.unrestricted(
                density_matrix_alph, density_matrix_beta, fock_matrix_alph, fock_matrix_beta
            )
            delta_energy = previous_total_energy - total_energy
            previous_total_energy = total_energy
            print('SCF ENERGY: ' + str(total_energy) + ' a.u.')

            if abs(delta_energy) < self.threshold:
                break

            energies_alph, coefficients_alph = self.linear_algebra.diagonalize(fock_matrix_alph)
            energies_beta, coefficients_beta = self.linear_algebra.diagonalize(fock_matrix_beta)

        return total_energy, energies_alph, energies_beta, coefficients_alph, coefficients_beta


class BlockedUnrestrictedSCF(SelfConsistentField):

    def __init__(self, linear_algebra, electrons, multiplicity, overlap, hamiltonian_matrix_factory):
        super().__init__(linear_algebra, electrons, hamiltonian_matrix_factory)
        self.calculate = TotalEnergy(hamiltonian_matrix_factory.core_hamiltonian)
        self.electrons_alph = (electrons + multiplicity - 1) // 2
        self.electrons_beta = (electrons - multiplicity + 1) // 2
        self.diis = DIIS(overlap, linear_algebra)

    def begin_iterations(self, orbital_energies, orbital_coefficients):
        total_energy = previous_total_energy = 0

        while True:

            if total_energy == 0:
                density_matrix = blocked_density_matrix(orbital_coefficients, orbital_coefficients.shape[0] // 2, 0)
            else:
                density_matrix = blocked_density_matrix(orbital_coefficients, self.electrons_alph, self.electrons_beta)

            fock_matrix = self.hamiltonian_matrix_factory.create(density_matrix)
            total_energy = self.calculate.restricted(density_matrix, fock_matrix)
            delta_energy = previous_total_energy - total_energy
            previous_total_energy = total_energy
            print('SCF ENERGY: ' + str(total_energy) + ' a.u.')

            if abs(delta_energy) < self.threshold:
                break

            fock_matrix = self.diis.fock_matrix(fock_matrix, density_matrix)
            orbital_energies, orbital_coefficients = self.linear_algebra.diagonalize(fock_matrix)

        return total_energy, orbital_energies, orbital_coefficients
