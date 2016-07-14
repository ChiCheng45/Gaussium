class TotalEnergy:
    """Calculates the total energy for hartree fock methods.

    Attributes
    ----------
    core_hamiltonian : np.matrix

    """

    def __init__(self, core_hamiltonian):
        self.core_hamiltonian = core_hamiltonian

    def restricted(self, density_matrix, fock_matrix):
        """Calculates the total energy for restricted hartree fock methods.

        Parameters
        ----------
        density_matrix : np.matrix
        fock_matrix : np.matrix

        Returns
        -------
        total_energy : float

        """
        length = density_matrix.shape[0]
        total_energy = 0
        for i in range(length):
            for j in range(length):
                total_energy += 1/2 * density_matrix.item(i, j) * (self.core_hamiltonian.item(i, j)
                + fock_matrix.item(i, j))
        return total_energy

    def unrestricted(self, density_matrix_alpha, density_matrix_beta, fock_matrix_alpha, fock_matrix_beta):
        """Calculates the total energy for unrestricted hartree fock methods.

        Parameters
        ----------
        density_matrix_alpha : np.matrix
        density_matrix_beta : np.matrix
        fock_matrix_alpha : np.matrix
        fock_matrix_beta : np.matrix

        Returns
        -------
        total_energy : float

        """
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
