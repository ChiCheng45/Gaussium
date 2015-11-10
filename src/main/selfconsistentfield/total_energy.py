class TotalEnergy:

    @staticmethod
    def calculate_total_energy(density_matrix, h_core_matrix, fock_matrix):
        a = density_matrix.shape[0]
        total_energy = 0
        for i in range(0, a):
            for j in range(0, a):
                total_energy += (1/2) * density_matrix.item((i, j)) * (h_core_matrix.item((i, j)) + fock_matrix.item((i, j)))
        return total_energy
