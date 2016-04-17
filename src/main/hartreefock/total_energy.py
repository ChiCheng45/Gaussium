import numpy as np


class TotalEnergy:

    def __init__(self, core_hamiltonian):
        self.core_hamiltonian = core_hamiltonian

    def restricted(self, density_matrix, fock_matrix):
        length = density_matrix.shape[0]
        total_energy = 0
        for i in range(length):
            for j in range(length):
                total_energy += (1/2) * density_matrix.item(i, j) * (self.core_hamiltonian.item(i, j)
                        + fock_matrix.item(i, j))
        return total_energy

    def unrestricted(self, density_matrix_alpha, density_matrix_beta, fock_matrix_alpha, fock_matrix_beta):
        density_matrix_total = density_matrix_alpha + density_matrix_beta
        length = density_matrix_total.shape[0]
        total_energy = 0
        for i in range(length):
            for j in range(length):
                out1 = density_matrix_total.item(i, j) * self.core_hamiltonian.item(i, j)
                out2 = density_matrix_alpha.item(i, j) * fock_matrix_alpha.item(i, j)
                out3 = density_matrix_beta.item(i, j) * fock_matrix_beta.item(i, j)
                total_energy += (1/2) * (out1 + out2 + out3)
        return total_energy

    def restricted_open_shell(self, density_matrix_doubly, density_matrix_singly, fock_matrix_doubly, fock_matrix_singly):
        out1 = np.trace((fock_matrix_doubly + self.core_hamiltonian) * density_matrix_doubly)
        out2 = 1/2 * np.trace((fock_matrix_singly + self.core_hamiltonian) * density_matrix_singly)
        total_energy = out1 + out2
        return total_energy
