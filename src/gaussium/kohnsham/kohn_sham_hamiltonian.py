import itertools
from gaussium.matrixelements import Matrix


class RestrictedKohnShamHamiltonian(Matrix):

    def __init__(self, core_hamiltonian, repulsion_matrix, exchange_correlation):
        super().__init__(repulsion_matrix.shape[0])
        self.core_hamiltonian = core_hamiltonian
        self.repulsion_matrix = repulsion_matrix
        self.exchange_correlation = exchange_correlation

    def create(self, density_matrix):

        def calculate_elst(i, j):
            g_ij = 0
            for a, b in itertools.product(range(self.matrix_size), repeat=2):
                g_ij += density_matrix.item(a, b) * self.repulsion_matrix.item(i, j, a, b)
            return g_ij

        def calculate_xc(i, j):
            return self.exchange_correlation.integrate_potential(density_matrix, i, j)

        return self.core_hamiltonian, self.create_matrix(calculate_elst), self.create_matrix(calculate_xc)
