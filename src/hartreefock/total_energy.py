import itertools


class TotalEnergy:
    """Calculates the total energy for hartree fock methods.

    Attributes
    ----------
    core_hamiltonian : np.matrix

    """

    def __init__(self, core_hamiltonian):
        self.core_hamiltonian = core_hamiltonian

    def restricted(self, density_matrix, hamiltonian):
        """Calculates the total energy for restricted hartree fock methods.

        Parameters
        ----------
        density_matrix : np.matrix
        hamiltonian : np.matrix

        Returns
        -------
        total_energy : float

        """
        length = density_matrix.shape[0]
        total_energy = 0
        for i, j in itertools.product(range(length), repeat=2):
            total_energy += 1/2 * density_matrix.item(i, j) * (self.core_hamiltonian.item(i, j)
            + hamiltonian.item(i, j))
        return total_energy

    def unrestricted(self, density_matrix_alpha, density_matrix_beta, hamiltonian_alpha, hamiltonian_beta):
        """Calculates the total energy for unrestricted hartree fock methods.

        Parameters
        ----------
        density_matrix_alpha : np.matrix
        density_matrix_beta : np.matrix
        hamiltonian_alpha : np.matrix
        hamiltonian_beta : np.matrix

        Returns
        -------
        total_energy : float

        """
        density_matrix_total = density_matrix_alpha + density_matrix_beta
        length = density_matrix_total.shape[0]
        total_energy = 0
        for i, j in itertools.product(range(length), repeat=2):
            out1 = density_matrix_total.item(i, j) * self.core_hamiltonian.item(i, j)
            out2 = density_matrix_alpha.item(i, j) * hamiltonian_alpha.item(i, j)
            out3 = density_matrix_beta.item(i, j) * hamiltonian_beta.item(i, j)
            total_energy += (1/2) * (out1 + out2 + out3)
        return total_energy
