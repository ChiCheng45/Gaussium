from unittest import TestCase
from src.main.matrixelements import DensityMatrixRestricted
from numpy import testing
import numpy as np


class TestDensityMatrixRestricted(TestCase):

    def setUp(self):
        orbital_coefficients = np.matrix([
            [0.9291, -0.6259],
            [0.1398,  1.1115]
        ])

        self.density_matrix = DensityMatrixRestricted(2)
        self.density_matrix.orbital_coefficient = orbital_coefficients
        self.density_matrix.matrix_size = 2

    def test_method_calculate_returns_density_matrix_element_for_0_0(self):
        element = self.density_matrix.calculate_restricted(0, 0)
        testing.assert_approx_equal(element, 1.7266, 4)

    def test_method_calculate_returns_density_matrix_element_for_1_0(self):
        element = self.density_matrix.calculate_restricted(1, 0)
        testing.assert_approx_equal(element, 0.2599, 3)

    def test_method_calculate_returns_density_matrix_element_for_1_1(self):
        element = self.density_matrix.calculate_restricted(1, 1)
        testing.assert_approx_equal(element, 0.0391, 3)
