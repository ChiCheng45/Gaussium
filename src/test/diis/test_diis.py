from src.main.diis_method import DIIS
from unittest import TestCase
from numpy import testing
import numpy as np


class TestBinomialCoefficientsFunction(TestCase):

    def setUp(self):
        orbital_coefficients = np.matrix([
            [1., 0.45077116],
            [0.45077116,  1.]
        ])

        self.diis = DIIS(orbital_coefficients)

    def test_fock_matrix_return_same_matrix_for_first_input(self):
        density_matrix = np.matrix([
            [1.72662514, 0.25985227],
            [0.25985227, 0.03910704]
        ])

        fock_matrix = np.matrix([
            [-1.3904162, -0.97320371],
            [-0.97320371, -0.74287664]
        ])

        diis_fock_matrix = self.diis.fock_matrix(fock_matrix, density_matrix)

        testing.assert_array_equal(diis_fock_matrix, fock_matrix)
