from src.main.common import Vector, BoysFunction
from math import sqrt, pi, exp
from src.main.objects import PrimitiveBasis


class ObaraSaika:

    def __init__(self):
        self.end_dict = {}

    def os_set(self, g1, g2, g3, g4):
        self.end_dict = {}
        return self.os_begin(g1, g2, g3, g4, 0)

    def os_begin(self, g1, g2, g3, g4, m):
        l_1 = g1.integral_exponents
        l_2 = g2.integral_exponents
        l_3 = g3.integral_exponents
        l_4 = g4.integral_exponents

        if l_1[0] + l_1[1] + l_1[2] + l_2[0] + l_2[1] + l_2[2] + l_3[0] + l_3[1] + l_3[2] + l_4[0] + l_4[1] + l_4[2] == 0:
            return self.os_end(g1, g2, g3, g4, m)
        else:
            if l_1[0] > 0:
                recursive_array = self.os_gaussian_factory(g1, g2, g3, g4, 0)
                return self.os_recursive(recursive_array, 0, m)
            if l_1[1] > 0:
                recursive_array = self.os_gaussian_factory(g1, g2, g3, g4, 1)
                return self.os_recursive(recursive_array, 1, m)
            if l_1[2] > 0:
                recursive_array = self.os_gaussian_factory(g1, g2, g3, g4, 2)
                return self.os_recursive(recursive_array, 2, m)
            if l_2[0] > 0:
                recursive_array = self.os_gaussian_factory(g2, g1, g4, g3, 0)
                return self.os_recursive(recursive_array, 0, m)
            if l_2[1] > 0:
                recursive_array = self.os_gaussian_factory(g2, g1, g4, g3, 1)
                return self.os_recursive(recursive_array, 1, m)
            if l_2[2] > 0:
                recursive_array = self.os_gaussian_factory(g2, g1, g4, g3, 2)
                return self.os_recursive(recursive_array, 2, m)
            if l_3[0] > 0:
                recursive_array = self.os_gaussian_factory(g3, g4, g1, g2, 0)
                return self.os_recursive(recursive_array, 0, m)
            if l_3[1] > 0:
                recursive_array = self.os_gaussian_factory(g3, g4, g1, g2, 1)
                return self.os_recursive(recursive_array, 1, m)
            if l_3[2] > 0:
                recursive_array = self.os_gaussian_factory(g3, g4, g1, g2, 2)
                return self.os_recursive(recursive_array, 2, m)
            if l_4[0] > 0:
                recursive_array = self.os_gaussian_factory(g4, g3, g2, g1, 0)
                return self.os_recursive(recursive_array, 0, m)
            if l_4[1] > 0:
                recursive_array = self.os_gaussian_factory(g4, g3, g2, g1, 1)
                return self.os_recursive(recursive_array, 1, m)
            if l_4[2] > 0:
                recursive_array = self.os_gaussian_factory(g4, g3, g2, g1, 2)
                return self.os_recursive(recursive_array, 2, m)

    def os_recursive(self, g, xyz, m):
        a_1 = g[0].exponent
        a_2 = g[1].exponent
        a_3 = g[2].exponent
        a_4 = g[3].exponent
        a_5 = a_1 + a_2
        a_6 = a_3 + a_4
        rho = (a_5 * a_6) / (a_5 + a_6)

        r_1 = g[0].coordinates
        r_2 = g[1].coordinates
        r_3 = g[2].coordinates
        r_4 = g[3].coordinates
        r_5 = Vector.gaussian(a_1, r_1, a_2, r_2)
        r_6 = Vector.gaussian(a_3, r_3, a_4, r_4)
        r_7 = Vector.gaussian(a_5, r_5, a_6, r_6)

        out1 = (r_5[xyz] - r_1[xyz]) * self.os_begin(g[0], g[1], g[2], g[3], m)
        out2 = (r_7[xyz] - r_5[xyz]) * self.os_begin(g[0], g[1], g[2], g[3], (m+1))
        if g[4].integral_exponents[xyz] >= 0:
            out3 = self.os_count(g[0].integral_exponents[xyz]) * (1 / (2 * a_5)) * self.os_begin(g[4], g[1], g[2], g[3], m)
            out4 = - self.os_count(g[0].integral_exponents[xyz]) * (rho / (2 * a_5**2)) * self.os_begin(g[4], g[1], g[2], g[3], (m+1))
        else:
            out3 = 0
            out4 = 0
        if g[5].integral_exponents[xyz] >= 0:
            out5 = self.os_count(g[1].integral_exponents[xyz]) * (1 / (2 * a_5)) * self.os_begin(g[0], g[5], g[2], g[3], m)
            out6 = - self.os_count(g[1].integral_exponents[xyz]) * (rho / (2 * a_5**2)) * self.os_begin(g[0], g[5], g[2], g[3], (m+1))
        else:
            out5 = 0
            out6 = 0
        if g[6].integral_exponents[xyz] >= 0:
            out7 = self.os_count(g[2].integral_exponents[xyz]) * (1 / (2*(a_5 + a_6))) * self.os_begin(g[0], g[1], g[6], g[3], (m+1))
        else:
            out7 = 0
        if g[7].integral_exponents[xyz] >= 0:
            out8 = self.os_count(g[3].integral_exponents[xyz]) * (1 / (2*(a_5 + a_6))) * self.os_begin(g[0], g[1], g[2], g[7], (m+1))
        else:
            out8 = 0
        return out1 + out2 + out3 + out4 + out5 + out6 + out7 + out8

    def os_end(self, g1, g2, g3, g4, m):
        if m in self.end_dict:
            return self.end_dict[m]
        else:
            a_1 = g1.exponent
            a_2 = g2.exponent
            a_3 = g3.exponent
            a_4 = g4.exponent
            a_5 = a_1 + a_2
            a_6 = a_3 + a_4

            r_1 = g1.coordinates
            r_2 = g2.coordinates
            r_3 = g3.coordinates
            r_4 = g4.coordinates
            r_5 = Vector.gaussian(a_1, r_1, a_2, r_2)
            r_6 = Vector.gaussian(a_3, r_3, a_4, r_4)
            r_12 = Vector.distance(r_1, r_2)
            r_34 = Vector.distance(r_3, r_4)
            r_56 = Vector.distance(r_5, r_6)

            t = (a_5 * a_6 * r_56**2) / (a_5 + a_6)
            out1 = (2 * pi**(5/2)) / (a_5 * a_6 * sqrt(a_5 + a_6))
            out2 = exp(((- a_1 * a_2 * r_12**2) / a_5) - ((a_3 * a_4 * r_34**2) / a_6))
            out3 = BoysFunction.function(m, t)
            ans = out1 * out2 * out3
            self.end_dict[m] = ans
        return ans

    @staticmethod
    def os_count(i):
        if i == 0:
            return 1
        else:
            return i

    @staticmethod
    def os_gaussian_factory(g1, g2, g3, g4, xyz):
        d_1 = g1.contraction
        d_2 = g2.contraction
        d_3 = g3.contraction
        d_4 = g4.contraction

        a_1 = g1.exponent
        a_2 = g2.exponent
        a_3 = g3.exponent
        a_4 = g4.exponent

        r_1 = g1.coordinates
        r_2 = g2.coordinates
        r_3 = g3.coordinates
        r_4 = g4.coordinates

        l_1 = g1.integral_exponents
        l_2 = g2.integral_exponents
        l_3 = g3.integral_exponents
        l_4 = g4.integral_exponents

        if xyz == 0:
            g1x1 = PrimitiveBasis(d_1, a_1, r_1, (l_1[0] - 1, l_1[1], l_1[2]))
            g1x2 = PrimitiveBasis(d_1, a_1, r_1, (l_1[0] - 2, l_1[1], l_1[2]))
            g2x1 = PrimitiveBasis(d_2, a_2, r_2, (l_2[0] - 1, l_2[1], l_2[2]))
            g3x1 = PrimitiveBasis(d_3, a_3, r_3, (l_3[0] - 1, l_3[1], l_3[2]))
            g4x1 = PrimitiveBasis(d_4, a_4, r_4, (l_4[0] - 1, l_4[1], l_4[2]))
            return g1x1, g2, g3, g4, g1x2, g2x1, g3x1, g4x1
        elif xyz == 1:
            g1y1 = PrimitiveBasis(d_1, a_1, r_1, (l_1[0], l_1[1] - 1, l_1[2]))
            g1y2 = PrimitiveBasis(d_1, a_1, r_1, (l_1[0], l_1[1] - 2, l_1[2]))
            g2y1 = PrimitiveBasis(d_2, a_2, r_2, (l_2[0], l_2[1] - 1, l_2[2]))
            g3y1 = PrimitiveBasis(d_3, a_3, r_3, (l_3[0], l_3[1] - 1, l_3[2]))
            g4y1 = PrimitiveBasis(d_4, a_4, r_4, (l_4[0], l_4[1] - 1, l_4[2]))
            return g1y1, g2, g3, g4, g1y2, g2y1, g3y1, g4y1
        elif xyz == 2:
            g1z1 = PrimitiveBasis(d_1, a_1, r_1, (l_1[0], l_1[1], l_1[2] - 1))
            g1z2 = PrimitiveBasis(d_1, a_1, r_1, (l_1[0], l_1[1], l_1[2] - 2))
            g2z1 = PrimitiveBasis(d_2, a_2, r_2, (l_2[0], l_2[1], l_2[2] - 1))
            g3z1 = PrimitiveBasis(d_3, a_3, r_3, (l_3[0], l_3[1], l_3[2] - 1))
            g4z1 = PrimitiveBasis(d_4, a_4, r_4, (l_4[0], l_4[1], l_4[2] - 1))
            return g1z1, g2, g3, g4, g1z2, g2z1, g3z1, g4z1
