from unittest import TestCase
from src.main.matrixelements import FockMatrixRestricted
from numpy import testing
import numpy as np


class TestGMatrixRestricted(TestCase):

    def setUp(self):
        repulsion_matrix = np.array([[[
            [1.3072, 0.4373],
            [0.4373, 0.6057]
        ], [
            [0.4373, 0.1773],
            [0.1773, 0.3118]
        ]], [[
            [0.4373, 0.1773],
            [0.1773, 0.3118]
        ], [
            [0.6057, 0.3118],
            [0.3118, 0.7746]
        ]]])

        density_matrix = np.matrix([
            [1.7266, 0.2599],
            [0.2599, 0.0391]
        ])

        self.g_matrix = FockMatrixRestricted(repulsion_matrix)
        self.g_matrix.density_matrix_total = density_matrix
        self.g_matrix.matrix_size = 2

    # def test_method_calculate_returns_g_matrix_element_for_0_0(self):
    #     element = self.g_matrix.calculate_restricted(0, 0)
    #     testing.assert_approx_equal(element, 1.2623, 4)
    #
    # def test_method_calculate_returns_g_matrix_element_for_1_0(self):
    #     element = self.g_matrix.calculate_restricted(1, 0)
    #     testing.assert_approx_equal(element, 0.3740, 4)
    #
    # def test_method_calculate_returns_g_matrix_element_for_1_1(self):
    #     element = self.g_matrix.calculate_restricted(1, 1)
    #     testing.assert_approx_equal(element, 0.9890, 4)
