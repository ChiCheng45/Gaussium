from src.main.matrixelements import Matrix
import itertools


class FockMatrixRestricted(Matrix):

    def __init__(self, core_hamiltonian, repulsion_matrix):
        super().__init__(repulsion_matrix.shape[0])
        self.repulsion_matrix = repulsion_matrix
        self.core_hamiltonian = core_hamiltonian

    def create(self, density_matrix):

        def calculate_restricted(i, j):
            g_ij = 0
            for a, b in itertools.product(range(self.matrix_size), repeat=2):
                coulomb_integral = self.repulsion_matrix.item(i, j, a, b)
                exchange_integral = self.repulsion_matrix.item(i, b, a, j)
                g_ij += density_matrix.item(a, b) * (coulomb_integral - 1/2 * exchange_integral)
            return g_ij

        return self.core_hamiltonian + self.create_matrix(calculate_restricted)


class FockMatrixUnrestricted(Matrix):

    def __init__(self, core_hamiltonian, repulsion_matrix):
        super().__init__(repulsion_matrix.shape[0])
        self.repulsion_matrix = repulsion_matrix
        self.core_hamiltonian = core_hamiltonian

    def create(self, density_matrix_alph, density_matrix_beta):
        density_matrix_total = density_matrix_alph + density_matrix_beta

        def unrestricted_alph(i, j):
            g_ij = 0
            for a, b in itertools.product(range(self.matrix_size), repeat=2):
                coulomb = density_matrix_total.item(a, b) * self.repulsion_matrix.item(i, j, a, b)
                exchange = density_matrix_alph.item(a, b) * self.repulsion_matrix.item(i, b, a, j)
                g_ij += coulomb - exchange
            return g_ij

        def unrestricted_beta(i, j):
            g_ij = 0
            for a, b in itertools.product(range(self.matrix_size), repeat=2):
                coulomb = density_matrix_total.item(a, b) * self.repulsion_matrix.item(i, j, a, b)
                exchange = density_matrix_beta.item(a, b) * self.repulsion_matrix.item(i, b, a, j)
                g_ij += coulomb - exchange
            return g_ij

        fock_matrix_alph = self.core_hamiltonian + self.create_matrix(unrestricted_alph)
        fock_matrix_beta = self.core_hamiltonian + self.create_matrix(unrestricted_beta)
        return fock_matrix_alph, fock_matrix_beta


class FockMatrixConstrained(FockMatrixUnrestricted):

    def __init__(self, core_hamiltonian, repulsion_matrix, electrons, multiplicity, linear_algebra):
        super().__init__(core_hamiltonian, repulsion_matrix)
        self.linear_algebra = linear_algebra
        self.electrons_alph = (electrons + multiplicity - 1) // 2
        self.electrons_beta = (electrons - multiplicity + 1) // 2

    def create(self, density_matrix_alph, density_matrix_beta):
        fock_matrix_alph, fock_matrix_beta = self.create(density_matrix_alph, density_matrix_beta)
        orthonormal_fock_alph = self.linear_algebra.orthonormalize(fock_matrix_alph)
        orthonormal_fock_beta = self.linear_algebra.orthonormalize(fock_matrix_beta)
        orthonormal_fock_closed = (orthonormal_fock_alph + orthonormal_fock_beta) / 2

        def constrained_alph(i, j):
            # (i is virtual and j is core) or (j is virtual and i is core)
            if i > (self.electrons_alph + 1 and self.electrons_beta + 1) > j \
            or j > (self.electrons_alph + 1 and self.electrons_beta + 1) > i:
                return orthonormal_fock_closed.item(i, j)
            else:
                return orthonormal_fock_alph.item(i, j)

        def constrained_beta(i, j):
            # (i is virtual and j is core) or (j is virtual and i is core)
            if i > (self.electrons_alph + 1 and self.electrons_beta + 1) > j \
            or j > (self.electrons_alph + 1 and self.electrons_beta + 1) > i:
                return orthonormal_fock_closed.item(i, j)
            else:
                return orthonormal_fock_beta.item(i, j)

        orthonormal_constrained_alph = self.create_matrix(constrained_alph)
        orthonormal_constrained_beta = self.create_matrix(constrained_beta)
        constrained_fock_matrix_alph = self.linear_algebra.non_orthonormal(orthonormal_constrained_alph)
        constrained_fock_matrix_beta = self.linear_algebra.non_orthonormal(orthonormal_constrained_beta)
        return constrained_fock_matrix_alph, constrained_fock_matrix_beta


class BlockedFockMatrixUnrestricted(Matrix):

    def __init__(self, core_hamiltonian, repulsion_matrix):
        super().__init__(repulsion_matrix.shape[0])
        self.repulsion_matrix = repulsion_matrix
        self.core_hamiltonian = core_hamiltonian

    def create(self, density_matrix):

        def coulomb(i, j):
            j_ij = 0
            for a, b in itertools.product(range(self.matrix_size), repeat=2):
                coulomb_integral = self.repulsion_matrix.item(i, j, a, b)
                j_ij += density_matrix.item(a, b) * coulomb_integral
            return j_ij

        def exchange(i, j):
            j_ij = 0
            for a, b in itertools.product(range(self.matrix_size), repeat=2):
                exchange_integral = self.repulsion_matrix.item(i, b, a, j)
                j_ij += density_matrix.item(a, b) * exchange_integral
            return j_ij

        coulomb = self.create_matrix(coulomb)
        exchange = self.create_matrix(exchange)
        return self.core_hamiltonian + coulomb - exchange
