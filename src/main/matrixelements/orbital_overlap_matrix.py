from src.main.integrals import OrbitalOverlap
from src.main.matrixelements import Matrix


class OrbitalOverlapMatrix(Matrix):

    def __init__(self, basis_set_array):
        super().__init__()
        self.basis_set_array = basis_set_array

    def create(self):
        self.matrix_size = len(self.basis_set_array)
        return self.create_matrix(self.calculate)

    def calculate(self, i, j):
        s_ij = 0
        if i == j:
            return 1
        else:
            primitive_array_i = self.basis_set_array[i].primitive_gaussian_array
            primitive_array_j = self.basis_set_array[j].primitive_gaussian_array
            for primitive_a in primitive_array_i:
                for primitive_b in primitive_array_j:
                    c_1 = primitive_a.contraction
                    c_2 = primitive_b.contraction
                    n_1 = primitive_a.normalisation
                    n_2 = primitive_b.normalisation
                    s_ij += n_1 * n_2 * c_1 * c_2 * OrbitalOverlap.integral(primitive_a, primitive_b)
            return s_ij
