from src.main.common import BinomialCoefficientsFunction
from scipy.misc import factorial2
from scipy.spatial import distance
import numpy as np


class OverlapIntegral:

    def s_function(self, l_1, l_2, a, b, gamma):
        s = 0
        for j in range(0, int((l_1 + l_2) / 2) + 1):
            s += BinomialCoefficientsFunction().calculate_coefficient(2*j, l_1, l_2, a, b) * (factorial2(2*j - 1) / (2*gamma)**j)
        return s

    def primitive_overlap_integral(self, gaussian_1, gaussian_2, r_a, r_b):
        r_ab = distance.euclidean(r_a, r_b)
        a_1 = gaussian_1.get_exponent()
        a_2 = gaussian_2.get_exponent()
        l_1 = gaussian_1.get_integral_exponents()
        l_2 = gaussian_2.get_integral_exponents()

        r_p = (a_1 * r_a + a_2 * r_b) / (a_1 + a_2)
        r_pa_xyz = r_p - r_a
        r_pb_xyz = r_p - r_b

        g = a_1 + a_2

        s_x = self.s_function(l_1[0], l_2[0], r_pa_xyz.item(0), r_pb_xyz.item(0), g)
        s_y = self.s_function(l_1[1], l_2[1], r_pa_xyz.item(1), r_pb_xyz.item(1), g)
        s_z = self.s_function(l_1[2], l_2[2], r_pa_xyz.item(2), r_pb_xyz.item(2), g)
        s_ij = (np.pi / g)**(3/2) * np.exp(- a_1 * a_2 * r_ab**2 / g) * s_x * s_y * s_z
        return s_ij
