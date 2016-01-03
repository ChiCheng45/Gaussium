from src.main.common import TotalEnergy
from math import floor, ceil


class SelfConsistentField:

    def __init__(self, core_hamiltonian, density_matrix_factory, g_matrix_factory, linear_algebra, electrons, multiplicity):
        self.core_hamiltonian = core_hamiltonian
        self.density_matrix_factory = density_matrix_factory
        self.g_matrix_factory = g_matrix_factory
        self.linear_algebra = linear_algebra
        self.electrons = electrons
        self.multiplicity = multiplicity
        self.total_energy = 0
        self.previous_total_energy = 0
        self.delta_energy = 1

    def restricted(self, orbital_coefficients):
        orbital_energies = []
        while abs(self.delta_energy) > 1e-9:
            density_matrix = self.density_matrix_factory.create_restricted(self.electrons, orbital_coefficients)
            g_matrix = self.g_matrix_factory.create_restricted(density_matrix)
            fock_matrix = self.core_hamiltonian + g_matrix
            orbital_energies, orbital_coefficients = self.linear_algebra.diagonalize(fock_matrix)
            self.total_energy = TotalEnergy.restricted(density_matrix, self.core_hamiltonian, fock_matrix)
            self.delta_energy = self.previous_total_energy - self.total_energy
            self.previous_total_energy = self.total_energy
            print('SCF ENERGY: ' + str(self.total_energy) + ' a.u.')
        return self.total_energy, orbital_energies, orbital_coefficients

    def unrestricted(self, orbital_coefficients):
        difference = floor(self.multiplicity / 2)
        electrons_alpha = ceil(self.electrons / 2) + difference
        electrons_beta = floor(self.electrons / 2) - difference
        coefficients_alpha = orbital_coefficients
        coefficients_beta = orbital_coefficients
        energies_alpha = []
        energies_beta = []

        while abs(self.delta_energy) > 1e-9:

            if self.total_energy == 0:
                density_matrix_alpha = self.density_matrix_factory.create_unrestricted((electrons_alpha + 1), coefficients_alpha)
                density_matrix_beta = self.density_matrix_factory.create_unrestricted((electrons_beta - 1), coefficients_beta)
            else:
                density_matrix_alpha = self.density_matrix_factory.create_unrestricted(electrons_alpha, coefficients_alpha)
                density_matrix_beta = self.density_matrix_factory.create_unrestricted(electrons_beta, coefficients_beta)

            g_matrix_alpha = self.g_matrix_factory.create_unrestricted(density_matrix_alpha, density_matrix_beta)
            g_matrix_beta = self.g_matrix_factory.create_unrestricted(density_matrix_beta, density_matrix_alpha)
            fock_matrix_alpha = self.core_hamiltonian + g_matrix_alpha
            fock_matrix_beta = self.core_hamiltonian + g_matrix_beta
            energies_alpha, coefficients_alpha = self.linear_algebra.diagonalize(fock_matrix_alpha)
            energies_beta, coefficients_beta = self.linear_algebra.diagonalize(fock_matrix_beta)

            self.total_energy = TotalEnergy.unrestricted(density_matrix_alpha, density_matrix_beta, self.core_hamiltonian, fock_matrix_alpha, fock_matrix_beta)
            self.delta_energy = self.previous_total_energy - self.total_energy
            self.previous_total_energy = self.total_energy
            print('SCF ENERGY: ' + str(self.total_energy) + ' a.u.')

        return self.total_energy, energies_alpha, energies_beta, coefficients_alpha, coefficients_beta
