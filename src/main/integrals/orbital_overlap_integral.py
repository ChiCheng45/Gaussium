from scipy.misc import factorial2
from src.main.common import BinomialCoefficientsFunction
from src.main.common.vector_manipulation import VectorManipulation
import numpy as np


class OverlapIntegral:

    def s_function(self, l_1, l_2, a, b, g):
        s = 0
        for j in range(0, int((l_1 + l_2) / 2) + 1):
            s += BinomialCoefficientsFunction.calculate_coefficient(2 * j, l_1, l_2, a, b) * (factorial2(2 * j - 1) / (2 * g)**j)
        return s

    def primitive_overlap_integral(self, gaussian_1, gaussian_2):


        a_1 = gaussian_1.exponent
        a_2 = gaussian_2.exponent
        l_1 = gaussian_1.integral_exponents
        l_2 = gaussian_2.integral_exponents

        r_a = gaussian_1.coordinates
        r_b = gaussian_2.coordinates
        r_ab = VectorManipulation.squared_distance(r_a, r_b)

        r_p = VectorManipulation.vector_gaussian(a_1, r_a, a_2, r_b)
        r_p_a = VectorManipulation.vector_minus(r_p, r_a)
        r_p_b = VectorManipulation.vector_minus(r_p, r_b)

        g = a_1 + a_2

        s_x = self.s_function(l_1[0], l_2[0], r_p_a[0], r_p_b[0], g)
        s_y = self.s_function(l_1[1], l_2[1], r_p_a[1], r_p_b[1], g)
        s_z = self.s_function(l_1[2], l_2[2], r_p_a[2], r_p_b[2], g)
        s_ij = (np.pi / g)**(3/2) * np.exp(- a_1 * a_2 * r_ab**2 / g) * s_x * s_y * s_z
        return s_ij
