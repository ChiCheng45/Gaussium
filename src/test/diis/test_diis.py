from unittest import TestCase
from src.main.diis import DIIS
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

    def test_fock_matrix_return_same_matrix_for_second_input(self):
        density_matrix_1 = np.matrix([
            [1.72662514, 0.25985227],
            [0.25985227, 0.03910704]
        ])

        fock_matrix_1 = np.matrix([
            [-1.3904162, -0.97320371],
            [-0.97320371, -0.74287664]
        ])

        density_matrix_2 = np.matrix([
            [1.33420227, 0.51662068],
            [0.51662068, 0.20004233]
        ])

        fock_matrix_2 = np.matrix([
            [-1.45139967, -1.04340121],
            [-1.04340121, -0.80339663]
        ])

        self.diis.fock_matrix(fock_matrix_1, density_matrix_1)
        diis_fock_matrix = self.diis.fock_matrix(fock_matrix_2, density_matrix_2)

        testing.assert_array_equal(diis_fock_matrix, fock_matrix_2)
