from src.diismethod import DIIS
from src.matrixelements import density_matrix_restricted
from src.hartreefock.scf_procedure import SelfConsistentField
from .kohn_sham_energy import KSEnergy


class KSRestrictedSCF(SelfConsistentField):

    def __init__(self, linear_algebra, electrons, overlap, exchange_correlation,
                 hamiltonian_matrix_factory):
        super().__init__(linear_algebra, electrons, hamiltonian_matrix_factory)
        self.diis = DIIS(overlap, linear_algebra)
        self.calculate = KSEnergy(self.electrons, exchange_correlation)

    def begin_iterations(self, orbital_energies, orbital_coefficients):
        previous_total_energy = 0

        while True:

            density_matrix = density_matrix_restricted(orbital_coefficients, self.electrons)
            core_matrix, elst_matrix, xc_matrix = self.hamiltonian_matrix_factory.create(density_matrix)
            ks_matrix = core_matrix + elst_matrix + xc_matrix
            total_energy = self.calculate.restricted(orbital_energies, density_matrix, elst_matrix, xc_matrix)
            delta_energy = previous_total_energy - total_energy
            previous_total_energy = total_energy
            print('SCF ENERGY: ' + str(total_energy) + ' a.u.')

            if abs(delta_energy) < self.threshold:
                break

            ks_matrix = self.diis.fock_matrix(ks_matrix, density_matrix)
            orbital_energies, orbital_coefficients = self.linear_algebra.diagonalize(ks_matrix)

        return total_energy, orbital_energies, orbital_coefficients