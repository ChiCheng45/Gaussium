import itertools
from src.factory import del_operator
from src.integrals import orbital_overlap
from src.matrixelements import Matrix


class KineticEnergyMatrix(Matrix):

    def __init__(self, basis_set_array):
        super().__init__(len(basis_set_array))
        self.basis_set_array = basis_set_array

    def create(self):
        return self.create_matrix(self.calculate)

    def calculate(self, i, j):
        t_ij = 0
        n_i = self.basis_set_array[i].normalisation
        n_j = self.basis_set_array[j].normalisation
        primitives_i = self.basis_set_array[i].primitive_gaussian_array
        primitives_j = self.basis_set_array[j].primitive_gaussian_array
        for primitive_a, primitive_b in itertools.product(primitives_i, primitives_j):
            c_1 = primitive_a.contraction
            c_2 = primitive_b.contraction
            n_1 = primitive_a.normalisation
            n_2 = primitive_b.normalisation
            primitives_k = del_operator(primitive_b)
            s_ij = 0
            for primitive_c in primitives_k:
                c_3 = primitive_c.contraction
                s_ij += n_1 * n_2 * c_1 * c_2 * c_3 * orbital_overlap(primitive_a, primitive_c)
            t_ij += s_ij
        return n_i * n_j * t_ij
