from src.main.objects import PrimitiveBasisFactory
from src.main.integrals import OverlapIntegral


class KineticEnergyElement:

    def __init__(self, basis_set_array):
        self.basis_set_array = basis_set_array

    def calculate(self, i, j):
        t_ij = 0
        primitive_gaussian_array_i = self.basis_set_array[i].primitive_gaussian_array
        primitive_gaussian_array_j = self.basis_set_array[j].primitive_gaussian_array
        for a in range(len(primitive_gaussian_array_i)):
            for b in range(len(primitive_gaussian_array_j)):
                c_1 = primitive_gaussian_array_i[a].contraction
                c_2 = primitive_gaussian_array_j[b].contraction
                n_1 = primitive_gaussian_array_i[a].normalisation()
                n_2 = primitive_gaussian_array_j[b].normalisation()
                primitive_gaussian_array_k = PrimitiveBasisFactory.del_operator(primitive_gaussian_array_j[b])
                s_ij = 0
                for c in range(len(primitive_gaussian_array_k)):
                    c_3 = primitive_gaussian_array_k[c].contraction
                    s_ij += n_1 * n_2 * c_1 * c_2 * c_3 * OverlapIntegral.primitive_overlap_integral(primitive_gaussian_array_i[a], primitive_gaussian_array_k[c])
                t_ij += s_ij
        return t_ij
