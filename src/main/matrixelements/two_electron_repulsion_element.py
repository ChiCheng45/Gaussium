from src.main.integrals.two_electron_repulsion_integral import ElectronRepulsionIntegral
from src.main.common import Symmetry
from multiprocessing import Pool


class TwoElectronRepulsionElement:

    def __init__(self, basis_set_array):
        self.basis_set_array = basis_set_array

    def calculate(self, i):
        if Symmetry.check_sym_2(self.basis_set_array[i[0]], self.basis_set_array[i[1]], self.basis_set_array[i[2]], self.basis_set_array[i[3]]):
            primitive_gaussian_array_i = self.basis_set_array[i[0]].primitive_gaussian_array
            primitive_gaussian_array_j = self.basis_set_array[i[1]].primitive_gaussian_array
            primitive_gaussian_array_k = self.basis_set_array[i[2]].primitive_gaussian_array
            primitive_gaussian_array_l = self.basis_set_array[i[3]].primitive_gaussian_array
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

    # single process method
    def store_integrals(self):
        repulsion_dict = {}
        for a in range(len(self.basis_set_array)):
            for b in range(len(self.basis_set_array)):
                for c in range(len(self.basis_set_array)):
                    for d in range(len(self.basis_set_array)):
                        if not (a > b or c > d or a > c or (a == c and b > d)):
                            index = (a, b, c, d)
                            repulsion_dict[index] = self.calculate(index)
        return repulsion_dict

    # multiprocess methods
    def keys_parallel(self):
        dict_key = []
        for a in range(len(self.basis_set_array)):
            for b in range(len(self.basis_set_array)):
                for c in range(len(self.basis_set_array)):
                    for d in range(len(self.basis_set_array)):
                        if not (a > b or c > d or a > c or (a == c and b > d)):
                            dict_key.append((a, b, c, d))
        return dict_key

    def store_parallel(self, processes):
        keys = self.keys_parallel()
        pool = Pool(processes)
        values = pool.map(self.calculate, keys)
        repulsion_dict = dict(zip(keys, values))
        return repulsion_dict
