from src.main.matrixelements import GElementUnrestricted, GElementRestricted
from src.main.matrixelements import DensityElementUnrestricted, DensityElementRestricted
from src.main.common import TotalEnergy
from math import floor, ceil


class SelfConsistentField:

    def __init__(self, core_hamiltonian, repulsion_dictionary, diagonalize, create_matrix, electrons, multiplicity):
        self.core_hamiltonian = core_hamiltonian
        self.repulsion = repulsion_dictionary
        self.diagonalize = diagonalize
        self.create_matrix = create_matrix
        self.electrons = electrons
        self.multiplicity = multiplicity
        self.total_energy = 0
        self.previous_total_energy = 0
        self.delta_energy = 1

    def restricted(self, orbital_coefficients):
        orbital_energies = []
        while abs(self.delta_energy) > 1e-9:
            density_matrix = self.create_matrix(DensityElementRestricted(orbital_coefficients, self.electrons))
            g_matrix = self.create_matrix(GElementRestricted(density_matrix, self.repulsion))
            fock_matrix = self.core_hamiltonian + g_matrix
            orbital_energies, orbital_coefficients = self.diagonalize(fock_matrix)
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
                density_matrix_alpha = self.create_matrix(DensityElementUnrestricted(coefficients_alpha, (electrons_alpha + 1)))
                density_matrix_beta = self.create_matrix(DensityElementUnrestricted(coefficients_beta, (electrons_beta - 1)))
            else:
                density_matrix_alpha = self.create_matrix(DensityElementUnrestricted(coefficients_alpha, electrons_alpha))
                density_matrix_beta = self.create_matrix(DensityElementUnrestricted(coefficients_beta, electrons_beta))

            density_matrix_total = density_matrix_alpha + density_matrix_beta
            g_matrix_alpha = self.create_matrix(GElementUnrestricted(density_matrix_total, density_matrix_alpha, self.repulsion))
            g_matrix_beta = self.create_matrix(GElementUnrestricted(density_matrix_total, density_matrix_beta, self.repulsion))
            fock_matrix_alpha = self.core_hamiltonian + g_matrix_alpha
            fock_matrix_beta = self.core_hamiltonian + g_matrix_beta
            energies_alpha, coefficients_alpha = self.diagonalize(fock_matrix_alpha)
            energies_beta, coefficients_beta = self.diagonalize(fock_matrix_beta)

            self.total_energy = TotalEnergy.unrestricted(density_matrix_total, density_matrix_alpha, density_matrix_beta, self.core_hamiltonian, fock_matrix_alpha, fock_matrix_beta)
            self.delta_energy = self.previous_total_energy - self.total_energy
            self.previous_total_energy = self.total_energy
            print('SCF ENERGY: ' + str(self.total_energy) + ' a.u.')

        return self.total_energy, energies_alpha, energies_beta, coefficients_alpha, coefficients_beta
