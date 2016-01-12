from src.main.common import Symmetry
from multiprocessing import Pool
import numpy as np


class MolecularIntegrals:

    def __init__(self, repulsion, orbital_coefficients):
        self.repulsion = repulsion
        self.orbital_coefficients = orbital_coefficients
        self.matrix_size = orbital_coefficients.shape[0]

    def calculate(self, i, j, k, l):
        repulsion_rstu = 0
        c = self.orbital_coefficients
        for r in range(self.matrix_size):
            repulsion_r = c.item(r, i) * self.repulsion[r, :, :, :]
            for s in range(self.matrix_size):
                repulsion_rs = c.item(s, j) * repulsion_r[s, :, :]
                for t in range(self.matrix_size):
                    repulsion_rst = c.item(t, k) * repulsion_rs[t, :]
                    for u in range(self.matrix_size):
                        repulsion_rstu += c.item(u, l) * repulsion_rst.item(u)
        return repulsion_rstu

    def store_parallel(self, processes):
        keys = []
        for a in range(self.matrix_size):
            for b in range(self.matrix_size):
                for c in range(self.matrix_size):
                    for d in range(self.matrix_size):
                        if not (a > b or c > d or a > c or (a == c and b > d)):
                            keys.append((a, b, c, d))

        pool = Pool(processes)
        values = pool.starmap(self.calculate, keys)
        repulsion_dict = dict(zip(keys, values))

        molecular_repulsion_matrix = np.zeros((self.matrix_size, self.matrix_size, self.matrix_size, self.matrix_size))
        for a in range(self.matrix_size):
            for b in range(self.matrix_size):
                for c in range(self.matrix_size):
                    for d in range(self.matrix_size):
                        molecular_repulsion_matrix.itemset((a, b, c, d), repulsion_dict[Symmetry.sort(a, b, c, d)])

        return molecular_repulsion_matrix
