from src.main.integrals.two_electron_repulsion_integral import ElectronRepulsionIntegral
from src.main.common import Symmetry


class TwoElectronRepulsionElement:

    def __init__(self, basis_set_array):
        self.basis_set_array = basis_set_array

    def calculate(self, i, j, k, l):
        if Symmetry.check_sym_2(self.basis_set_array[i], self.basis_set_array[j], self.basis_set_array[k], self.basis_set_array[l]):
            primitive_gaussian_array_i = self.basis_set_array[i].primitive_gaussian_array
            primitive_gaussian_array_j = self.basis_set_array[j].primitive_gaussian_array
            primitive_gaussian_array_k = self.basis_set_array[k].primitive_gaussian_array
            primitive_gaussian_array_l = self.basis_set_array[l].primitive_gaussian_array
            f_mn = 0
            for a in range(len(primitive_gaussian_array_i)):
                for b in range(len(primitive_gaussian_array_j)):
                    for c in range(len(primitive_gaussian_array_k)):
                        for d in range(len(primitive_gaussian_array_l)):
                            c_1 = primitive_gaussian_array_i[a].contraction
                            c_2 = primitive_gaussian_array_j[b].contraction
                            c_3 = primitive_gaussian_array_k[c].contraction
                            c_4 = primitive_gaussian_array_l[d].contraction
                            n_1 = primitive_gaussian_array_i[a].normalisation
                            n_2 = primitive_gaussian_array_j[b].normalisation
                            n_3 = primitive_gaussian_array_k[c].normalisation
                            n_4 = primitive_gaussian_array_l[d].normalisation
                            integral = ElectronRepulsionIntegral.integral(primitive_gaussian_array_i[a], primitive_gaussian_array_j[b], primitive_gaussian_array_k[c], primitive_gaussian_array_l[d])
                            f_mn += c_1 * c_2 * c_3 * c_4 * n_1 * n_2 * n_3 * n_4 * integral
            return f_mn
        else:
            return 0

    def store_integrals(self):
        electron_repulsion_dict = {}
        for a in range(len(self.basis_set_array)):
            for b in range(len(self.basis_set_array)):
                for c in range(len(self.basis_set_array)):
                    for d in range(len(self.basis_set_array)):
                        if not (a > b or c > d or a > c or (a == c and b > d)):
                            electron_repulsion_dict[(a, b, c, d)] = self.calculate(a, b, c, d)
        return electron_repulsion_dict
