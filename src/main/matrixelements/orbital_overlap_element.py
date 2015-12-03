from src.main.integrals import OverlapIntegral


class OrbitalOverlapElement:

    def __init__(self, basis_set_array):
        self.basis_set_array = basis_set_array

    def calculate(self, i, j):
        s_ij = 0
        if i == j:
            return 1
        else:
            primitive_gaussian_array_i = self.basis_set_array[i].primitive_gaussian_array
            primitive_gaussian_array_j = self.basis_set_array[j].primitive_gaussian_array
            for a in range(len(primitive_gaussian_array_i)):
                for b in range(len(primitive_gaussian_array_j)):
                    c_1 = primitive_gaussian_array_i[a].contraction
                    c_2 = primitive_gaussian_array_j[b].contraction
                    n_1 = primitive_gaussian_array_i[a].normalisation()
                    n_2 = primitive_gaussian_array_j[b].normalisation()
                    s_ij += n_1 * n_2 * c_1 * c_2 * OverlapIntegral.primitive_overlap_integral(primitive_gaussian_array_i[a], primitive_gaussian_array_j[b])
            return s_ij
