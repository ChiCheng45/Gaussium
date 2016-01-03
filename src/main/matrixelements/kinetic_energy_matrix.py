from src.main.objects import PrimitiveBasisFactory
from src.main.integrals import OrbitalOverlap
from src.main.matrixelements import Matrix


class KineticEnergyMatrix(Matrix):

    def __init__(self):
        super().__init__()
        self.basis_set_array = []

    def create(self, basis_set_array):
        self.basis_set_array = basis_set_array
        self.matrix_size = len(basis_set_array)
        return self.create_matrix(self.calculate)

    def calculate(self, i, j):
        t_ij = 0
        primitive_array_i = self.basis_set_array[i].primitive_gaussian_array
        primitive_array_j = self.basis_set_array[j].primitive_gaussian_array
        for primitive_a in primitive_array_i:
            for primitive_b in primitive_array_j:
                c_1 = primitive_a.contraction
                c_2 = primitive_b.contraction
                n_1 = primitive_a.normalisation
                n_2 = primitive_b.normalisation
                primitive_gaussian_array_k = PrimitiveBasisFactory.del_operator(primitive_b)
                s_ij = 0
                for primitive_c in primitive_gaussian_array_k:
                    c_3 = primitive_c.contraction
                    s_ij += n_1 * n_2 * c_1 * c_2 * c_3 * OrbitalOverlap.integral(primitive_a, primitive_c)
                t_ij += s_ij
        return t_ij
