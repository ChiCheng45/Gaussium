from src.main.common import Vector
from src.main.integrals import BoysFunction
from math import sqrt, pi, exp
from src.main.objects import PrimitiveBasis


class ObaraSaika:

    def __init__(self):
        self.a_7 = 0
        self.r_7 = ()
        self.boys_out1 = 0
        self.boys_x = 0
        self.end_dict = {}

    def os_set(self, g1, g2, g3, g4):
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
        r_5 = Vector.gaussian_product(a_1, r_1, a_2, r_2)
        r_6 = Vector.gaussian_product(a_3, r_3, a_4, r_4)
        r_12 = Vector.distance(r_1, r_2)
        r_34 = Vector.distance(r_3, r_4)
        r_56 = Vector.distance(r_5, r_6)

        self.a_7 = (a_5 * a_6) / (a_5 + a_6)
        self.r_7 = Vector.gaussian_product(a_5, r_5, a_6, r_6)
        self.boys_out1 = (2 * pi**(5/2)) / (a_5 * a_6 * sqrt(a_5 + a_6)) * exp(((- a_1 * a_2 * r_12**2) / a_5) - ((a_3 * a_4 * r_34**2) / a_6))
        self.boys_x = (a_5 * a_6 * r_56**2) / (a_5 + a_6)
        self.end_dict = {}

        return self.os_begin(0, g1, g2, g3, g4)

    def os_begin(self, m, g1, g2, g3, g4):
        l_1 = g1.integral_exponents
        l_2 = g2.integral_exponents
        l_3 = g3.integral_exponents
        l_4 = g4.integral_exponents

        if l_1[0] > 0:
            return self.os_recursive(0, m, *self.os_gaussian_factory(0, g1, g2, g3, g4))
        elif l_1[1] > 0:
            return self.os_recursive(1, m, *self.os_gaussian_factory(1, g1, g2, g3, g4))
        elif l_1[2] > 0:
            return self.os_recursive(2, m, *self.os_gaussian_factory(2, g1, g2, g3, g4))
        elif l_2[0] > 0:
            return self.os_recursive(0, m, *self.os_gaussian_factory(0, g2, g1, g4, g3))
        elif l_2[1] > 0:
            return self.os_recursive(1, m, *self.os_gaussian_factory(1, g2, g1, g4, g3))
        elif l_2[2] > 0:
            return self.os_recursive(2, m, *self.os_gaussian_factory(2, g2, g1, g4, g3))
        elif l_3[0] > 0:
            return self.os_recursive(0, m, *self.os_gaussian_factory(0, g3, g4, g1, g2))
        elif l_3[1] > 0:
            return self.os_recursive(1, m, *self.os_gaussian_factory(1, g3, g4, g1, g2))
        elif l_3[2] > 0:
            return self.os_recursive(2, m, *self.os_gaussian_factory(2, g3, g4, g1, g2))
        elif l_4[0] > 0:
            return self.os_recursive(0, m, *self.os_gaussian_factory(0, g4, g3, g2, g1))
        elif l_4[1] > 0:
            return self.os_recursive(1, m, *self.os_gaussian_factory(1, g4, g3, g2, g1))
        elif l_4[2] > 0:
            return self.os_recursive(2, m, *self.os_gaussian_factory(2, g4, g3, g2, g1))
        else:
            return self.os_end(m)

    def os_recursive(self, r, m, g1, g2, g3, g4, g5, g6, g7, g8):
        out3 = out4 = out5 = out6 = out7 = out8 = 0

        a_1 = g1.exponent
        a_2 = g2.exponent
        a_3 = g3.exponent
        a_4 = g4.exponent
        a_5 = a_1 + a_2
        a_6 = a_3 + a_4

        r_1 = g1.coordinates
        r_2 = g2.coordinates
        r_5 = Vector.gaussian_product(a_1, r_1, a_2, r_2)

        out1 = (r_5[r] - r_1[r]) * self.os_begin(m, g1, g2, g3, g4)
        out2 = (self.r_7[r] - r_5[r]) * self.os_begin((m+1), g1, g2, g3, g4)
        if g5.integral_exponents[r] >= 0:
            out3 = self.os_count(g1.integral_exponents[r]) * (1 / (2 * a_5)) * self.os_begin(m, g5, g2, g3, g4)
            out4 = self.os_count(g1.integral_exponents[r]) * (self.a_7 / (2 * a_5**2)) * self.os_begin((m+1), g5, g2, g3, g4)
        if g6.integral_exponents[r] >= 0:
            out5 = self.os_count(g2.integral_exponents[r]) * (1 / (2 * a_5)) * self.os_begin(m, g1, g6, g3, g4)
            out6 = self.os_count(g2.integral_exponents[r]) * (self.a_7 / (2 * a_5**2)) * self.os_begin((m+1), g1, g6, g3, g4)
        if g7.integral_exponents[r] >= 0:
            out7 = self.os_count(g3.integral_exponents[r]) * (1 / (2*(a_5 + a_6))) * self.os_begin((m+1), g1, g2, g7, g4)
        if g8.integral_exponents[r] >= 0:
            out8 = self.os_count(g4.integral_exponents[r]) * (1 / (2*(a_5 + a_6))) * self.os_begin((m+1), g1, g2, g3, g8)

        return out1 + out2 + out3 - out4 + out5 - out6 + out7 + out8

    def os_count(self, i):
        if i == 0:
            return 1
        else:
            return i

    def os_gaussian_factory(self, r, g1, g2, g3, g4):
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

        if r == 0:
            g1x1 = PrimitiveBasis(d_1, a_1, r_1, (l_1[0] - 1, l_1[1], l_1[2]))
            g1x2 = PrimitiveBasis(d_1, a_1, r_1, (l_1[0] - 2, l_1[1], l_1[2]))
            g2x1 = PrimitiveBasis(d_2, a_2, r_2, (l_2[0] - 1, l_2[1], l_2[2]))
            g3x1 = PrimitiveBasis(d_3, a_3, r_3, (l_3[0] - 1, l_3[1], l_3[2]))
            g4x1 = PrimitiveBasis(d_4, a_4, r_4, (l_4[0] - 1, l_4[1], l_4[2]))
            return g1x1, g2, g3, g4, g1x2, g2x1, g3x1, g4x1
        elif r == 1:
            g1y1 = PrimitiveBasis(d_1, a_1, r_1, (l_1[0], l_1[1] - 1, l_1[2]))
            g1y2 = PrimitiveBasis(d_1, a_1, r_1, (l_1[0], l_1[1] - 2, l_1[2]))
            g2y1 = PrimitiveBasis(d_2, a_2, r_2, (l_2[0], l_2[1] - 1, l_2[2]))
            g3y1 = PrimitiveBasis(d_3, a_3, r_3, (l_3[0], l_3[1] - 1, l_3[2]))
            g4y1 = PrimitiveBasis(d_4, a_4, r_4, (l_4[0], l_4[1] - 1, l_4[2]))
            return g1y1, g2, g3, g4, g1y2, g2y1, g3y1, g4y1
        elif r == 2:
            g1z1 = PrimitiveBasis(d_1, a_1, r_1, (l_1[0], l_1[1], l_1[2] - 1))
            g1z2 = PrimitiveBasis(d_1, a_1, r_1, (l_1[0], l_1[1], l_1[2] - 2))
            g2z1 = PrimitiveBasis(d_2, a_2, r_2, (l_2[0], l_2[1], l_2[2] - 1))
            g3z1 = PrimitiveBasis(d_3, a_3, r_3, (l_3[0], l_3[1], l_3[2] - 1))
            g4z1 = PrimitiveBasis(d_4, a_4, r_4, (l_4[0], l_4[1], l_4[2] - 1))
            return g1z1, g2, g3, g4, g1z2, g2z1, g3z1, g4z1

    def os_end(self, m):
        if m in self.end_dict:
            return self.end_dict[m]
        else:
            boys_out2 = BoysFunction.calculate(m, self.boys_x)
            ans = self.boys_out1 * boys_out2
            self.end_dict[m] = ans
        return ans
