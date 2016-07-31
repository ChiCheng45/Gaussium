from src.main.integrals import ElectronRepulsion
from src.main.integrals import ObaraSaika
from src.main.integrals import HeadGordonPople
from multiprocessing import Pool
import numpy as np
import itertools


class TwoElectronRepulsion:

    def __init__(self, basis_set_array, integral, symmetry, processes):
        self.basis_set_array = basis_set_array
        self.matrix_size = len(basis_set_array)
        self.integral = integral
        self.symmetry = symmetry
        self.processes = processes

    def calculate(self, i, j, k, l):
        if self.symmetry.none_zero_integral((i, j, k, l)):
            f_mn = 0.0
            primitives_i = self.basis_set_array[i].primitive_gaussian_array
            primitives_j = self.basis_set_array[j].primitive_gaussian_array
            primitives_k = self.basis_set_array[k].primitive_gaussian_array
            primitives_l = self.basis_set_array[l].primitive_gaussian_array
            for primitive_a, primitive_b, primitive_c, primitive_d in itertools.product(primitives_i, primitives_j,
            primitives_k, primitives_l):
                c_1 = primitive_a.contraction
                c_2 = primitive_b.contraction
                c_3 = primitive_c.contraction
                c_4 = primitive_d.contraction
                n_1 = primitive_a.normalisation
                n_2 = primitive_b.normalisation
                n_3 = primitive_c.normalisation
                n_4 = primitive_d.normalisation
                integral = self.integral.integrate(primitive_a, primitive_b, primitive_c, primitive_d)
                f_mn += c_1 * c_2 * c_3 * c_4 * n_1 * n_2 * n_3 * n_4 * integral
            return f_mn
        else:
            return 0.0

    def create(self):

        keys = []
        for a, b, c, d in itertools.product(range(self.matrix_size), repeat=4):
            if not (a > b or c > d or a > c or (a == c and b > d)):
                keys.append((a, b, c, d))

        if self.processes > 1:
            pool = Pool(self.processes)
            values = pool.starmap(self.calculate, keys)
            pool.close()
            repulsion_dictionary = dict(zip(keys, values))
        else:
            repulsion_dictionary = {index: self.calculate(*index) for index in keys}

        repulsion_matrix = np.zeros((self.matrix_size, self.matrix_size, self.matrix_size, self.matrix_size))
        for a, b, c, d in itertools.product(range(self.matrix_size), repeat=4):
            repulsion_matrix.itemset((a, b, c, d), repulsion_dictionary[self.symmetry.sort_index(a, b, c, d)])

        return repulsion_matrix


class TwoElectronRepulsionMatrixOS(TwoElectronRepulsion):

    def __init__(self, basis_set_array, symmetry_matrix, processes):
        super().__init__(basis_set_array, ObaraSaika(), symmetry_matrix, processes)


class TwoElectronRepulsionMatrixCook(TwoElectronRepulsion):

    def __init__(self, basis_set_array, symmetry_matrix, processes):
        super().__init__(basis_set_array, ElectronRepulsion(), symmetry_matrix, processes)


class TwoElectronRepulsionMatrixHGP(TwoElectronRepulsion):

    def __init__(self, basis_set_array, symmetry_matrix, processes):
        super().__init__(basis_set_array, HeadGordonPople(), symmetry_matrix, processes)
