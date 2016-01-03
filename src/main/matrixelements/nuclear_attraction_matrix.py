from src.main.integrals import NuclearAttraction
from src.main.matrixelements import Matrix


class NuclearAttractionMatrix(Matrix):

    def __init__(self):
        super().__init__()
        self.nuclei_array = []
        self.basis_set_array = []

    def create(self, basis_set_array, nuclei_array):
        self.basis_set_array = basis_set_array
        self.nuclei_array = nuclei_array
        self.matrix_size = len(basis_set_array)
        return self.create_matrix(self.calculate)

    def calculate(self, i, j):
        v_ij = 0
        primitive_array_i = self.basis_set_array[i].primitive_gaussian_array
        primitive_array_j = self.basis_set_array[j].primitive_gaussian_array
        for primitive_a in primitive_array_i:
            for primitive_b in primitive_array_j:
                c_1 = primitive_a.contraction
                c_2 = primitive_b.contraction
                n_1 = primitive_a.normalisation
                n_2 = primitive_b.normalisation
                for nuclei in self.nuclei_array:
                    v_ij += - nuclei.charge * n_1 * n_2 * c_1 * c_2 * NuclearAttraction.integral(primitive_a, primitive_b, nuclei)
        return v_ij
