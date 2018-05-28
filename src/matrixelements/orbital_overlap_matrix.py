import itertools
from src.integrals import orbital_overlap
from src.matrixelements import Matrix


class OrbitalOverlapMatrix(Matrix):

    def __init__(self, basis_set_array):
        super().__init__(len(basis_set_array))
        self.basis_set_array = basis_set_array

    def create(self):
        return self.create_matrix(self.calculate)

    def calculate(self, i, j):
        s_ij = 0
        n_i = self.basis_set_array[i].normalisation
        n_j = self.basis_set_array[j].normalisation
        primitives_i = self.basis_set_array[i].primitive_gaussian_array
        primitives_j = self.basis_set_array[j].primitive_gaussian_array
        for primitive_a, primitive_b in itertools.product(primitives_i, primitives_j):
            c_1 = primitive_a.contraction
            c_2 = primitive_b.contraction
            n_1 = primitive_a.normalisation
            n_2 = primitive_b.normalisation
            s_ij += n_1 * n_2 * c_1 * c_2 * orbital_overlap(primitive_a, primitive_b)
        return n_i * n_j * s_ij
