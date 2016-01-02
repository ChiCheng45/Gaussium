from unittest import TestCase
from src.main.matrixelements import DensityElementRestricted
from numpy import testing
import numpy as np


class TestDensityElementRestricted(TestCase):

    def setUp(self):
        self.electrons = 2
        self.orbital_coefficient_matrix = np.matrix([[0.9291, -0.6259], [0.1398, 1.1115]])

    def test_method_calculate_returns_density_matrix_element_for_0_0(self):
        density_matrix_element = DensityElementRestricted(self.orbital_coefficient_matrix, self.electrons)
        element = density_matrix_element.calculate(0, 0)
        testing.assert_approx_equal(element, 1.7266, 4)

    def test_method_calculate_returns_density_matrix_element_for_1_0(self):
        density_matrix_element = DensityElementRestricted(self.orbital_coefficient_matrix, self.electrons)
        element = density_matrix_element.calculate(1, 0)
        testing.assert_approx_equal(element, 0.2599, 3)

    def test_method_calculate_returns_density_matrix_element_for_1_1(self):
        density_matrix_element = DensityElementRestricted(self.orbital_coefficient_matrix, self.electrons)
        element = density_matrix_element.calculate(1, 1)
        testing.assert_approx_equal(element, 0.0391, 3)
