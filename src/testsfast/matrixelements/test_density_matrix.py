from unittest import TestCase
from src.main.matrixelements import density_matrix_restricted
from numpy import testing
import numpy as np


class TestDensityMatrixRestricted(TestCase):

    def setUp(self):
        self.orbital_coefficients = np.matrix([
                [0.9291, -0.6259],
                [0.1398,  1.1115]
        ])

    def test_method_calculate_returns_density_matrix_element_for_0_0(self):
        element = density_matrix_restricted(self.orbital_coefficients, 2).item(0, 0)
        testing.assert_approx_equal(element, 1.7266, 4)

    def test_method_calculate_returns_density_matrix_element_for_1_0(self):
        element = density_matrix_restricted(self.orbital_coefficients, 2).item(1, 0)
        testing.assert_approx_equal(element, 0.2599, 3)

    def test_method_calculate_returns_density_matrix_element_for_1_1(self):
        element = density_matrix_restricted(self.orbital_coefficients, 2).item(1, 1)
        testing.assert_approx_equal(element, 0.0391, 3)
