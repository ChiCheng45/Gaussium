from src.main.matrixelements import Matrix
import itertools


class RestrictedKohnShamHamiltonian(Matrix):

    def __init__(self, core_hamiltonian, repulsion_matrix, exchange_integral):
        super().__init__(repulsion_matrix.shape[0])
        self.core_hamiltonian = core_hamiltonian
        self.repulsion_matrix = repulsion_matrix
        self.exchange_integral = exchange_integral

    def create(self, density_matrix):

        def calculate_restricted(i, j):
            g_ij = 0
            for a, b in itertools.product(range(self.matrix_size), repeat=2):
                coulomb = self.repulsion_matrix.item(i, j, a, b)
                g_ij += density_matrix.item(a, b) * coulomb
            g_ij += self.exchange_integral.integrate(density_matrix, i, j)
            return g_ij

        return self.core_hamiltonian + self.create_matrix(calculate_restricted)
