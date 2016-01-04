from multiprocessing import Pool
from src.main.integrals.twoelectronrepulsion import ElectronRepulsion, ObaraSaika, HeadGordonPople
from src.main.common import Symmetry


class TwoElectronRepulsionElement:

    def __init__(self, basis_set_array, integral):
        self.basis_set_array = basis_set_array
        self.integral = integral

    def calculate(self, i, j, k, l):
        if Symmetry.check_sym(self.basis_set_array[i], self.basis_set_array[j], self.basis_set_array[k], self.basis_set_array[l]):
            f_mn = 0
            primitive_array_i = self.basis_set_array[i].primitive_gaussian_array
            primitive_array_j = self.basis_set_array[j].primitive_gaussian_array
            primitive_array_k = self.basis_set_array[k].primitive_gaussian_array
            primitive_array_l = self.basis_set_array[l].primitive_gaussian_array
            for primitive_a in primitive_array_i:
                for primitive_b in primitive_array_j:
                    for primitive_c in primitive_array_k:
                        for primitive_d in primitive_array_l:
                            c_1 = primitive_a.contraction
                            c_2 = primitive_b.contraction
                            c_3 = primitive_c.contraction
                            c_4 = primitive_d.contraction
                            n_1 = primitive_a.normalisation
                            n_2 = primitive_b.normalisation
                            n_3 = primitive_c.normalisation
                            n_4 = primitive_d.normalisation
                            integral = self.integral(primitive_a, primitive_b, primitive_c, primitive_d)
                            f_mn += c_1 * c_2 * c_3 * c_4 * n_1 * n_2 * n_3 * n_4 * integral
            return f_mn
        else:
            return 0

    def store_series(self):
        repulsion_dict = {}
        for a in range(len(self.basis_set_array)):
            for b in range(len(self.basis_set_array)):
                for c in range(len(self.basis_set_array)):
                    for d in range(len(self.basis_set_array)):
                        if not (a > b or c > d or a > c or (a == c and b > d)):
                            index = (a, b, c, d)
                            repulsion_dict[index] = self.calculate(*index)
        return repulsion_dict

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
        if processes == 1:
            return self.store_series()
        else:
            keys = self.keys_parallel()
            pool = Pool(processes)
            values = pool.starmap(self.calculate, keys)
            repulsion_dict = dict(zip(keys, values))
            return repulsion_dict


class TwoElectronRepulsionElementCook(TwoElectronRepulsionElement):

    def __init__(self, basis_set_array):
        super().__init__(basis_set_array, ElectronRepulsion().integral)


class TwoElectronRepulsionElementOS(TwoElectronRepulsionElement):

    def __init__(self, basis_set_array):
        super().__init__(basis_set_array, ObaraSaika().os_set)


class TwoElectronRepulsionElementHGP(TwoElectronRepulsionElement):

    def __init__(self, basis_set_array):
        super().__init__(basis_set_array, HeadGordonPople().hgp_set)
