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

    def calculate_integral(self, i, j, k, l):
        if self.symmetry.none_zero_integral((i, j, k, l)):
            basis_i = self.basis_set_array[i]
            basis_j = self.basis_set_array[j]
            basis_k = self.basis_set_array[k]
            basis_l = self.basis_set_array[l]
            return self.integral.integrate(basis_i, basis_j, basis_k, basis_l)
        else:
            return 0.0

    def create_repulsion_matrix(self):

        keys = []
        for a, b, c, d in itertools.product(range(self.matrix_size), repeat=4):
            if not (a > b or c > d or a > c or (a == c and b > d)):
                keys.append((a, b, c, d))

        if self.processes > 1:
            pool = Pool(self.processes)
            values = pool.starmap(self.calculate_integral, keys)
            pool.close()
            repulsion_dictionary = dict(zip(keys, values))
        else:
            repulsion_dictionary = {index: self.calculate_integral(*index) for index in keys}

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
