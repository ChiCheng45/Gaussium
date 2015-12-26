class TotalEnergy:

    @staticmethod
    def restricted(density_matrix, core_hamiltonian, fock_matrix):
        a = density_matrix.shape[0]
        total_energy = 0
        for i in range(a):
            for j in range(a):
                total_energy += (1/2) * density_matrix.item(i, j) * (core_hamiltonian.item(i, j) + fock_matrix.item(i, j))
        return total_energy

    @staticmethod
    def unrestricted(density_matrix_total, density_matrix_alpha, density_matrix_beta, core_hamiltonian, fock_matrix_alpha, fock_matrix_beta):
        a = density_matrix_total.shape[0]
        total_energy = 0
        for i in range(a):
            for j in range(a):
                out1 = density_matrix_total.item(i, j) * core_hamiltonian.item(i, j)
                out2 = density_matrix_alpha.item(i, j) * fock_matrix_alpha.item(i, j)
                out3 = density_matrix_beta.item(i, j) * fock_matrix_beta.item(i, j)
                total_energy += (1/2) * (out1 + out2 + out3)
        return total_energy
