from src.main.common import Vector, BoysFunction
from math import sqrt, pi, exp
from src.main.objects import PrimitiveBasis


class HeadGordonPople:

    def __init__(self):
        self.end_dict = {}

    def hgp_set(self, g1, g2, g3, g4):
        self.end_dict = {}
        return self.hgp_sort(g1, g2, g3, g4)

    def hgp_sort(self, g1, g2, g3, g4):
        l_1 = g1.integral_exponents
        l_2 = g2.integral_exponents
        l_3 = g3.integral_exponents
        l_4 = g4.integral_exponents
        if sum(l_1) >= sum(l_2):
            if sum(l_3) >= sum(l_4):
                return self.hgp_begin_horizontal(g1, g2, g3, g4)
            else:
                return self.hgp_begin_horizontal(g1, g2, g4, g3)
        else:
            if sum(l_3) >= sum(l_4):
                return self.hgp_begin_horizontal(g2, g1, g3, g4)
            else:
                return self.hgp_begin_horizontal(g2, g1, g4, g3)

    def hgp_begin_horizontal(self, g1, g2, g3, g4):
        l_2 = g2.integral_exponents
        l_4 = g4.integral_exponents
        if l_2[0] > 0:
            gaussian_tuple = self.hgp_horizontal_factory(0, g1, g2, g3, g4)
            return self.horizontal_recursion(0, *gaussian_tuple)
        elif l_2[1] > 0:
            gaussian_tuple = self.hgp_horizontal_factory(1, g1, g2, g3, g4)
            return self.horizontal_recursion(1, *gaussian_tuple)
        elif l_2[2] > 0:
            gaussian_tuple = self.hgp_horizontal_factory(2, g1, g2, g3, g4)
            return self.horizontal_recursion(2, *gaussian_tuple)
        elif l_4[0] > 0:
            gaussian_tuple = self.hgp_horizontal_factory(0, g3, g4, g1, g2)
            return self.horizontal_recursion(0, *gaussian_tuple)
        elif l_4[1] > 0:
            gaussian_tuple = self.hgp_horizontal_factory(1, g3, g4, g1, g2)
            return self.horizontal_recursion(1, *gaussian_tuple)
        elif l_4[2] > 0:
                gaussian_tuple = self.hgp_horizontal_factory(2, g3, g4, g1, g2)
                return self.horizontal_recursion(2, *gaussian_tuple)
        else:
            return self.hgp_begin_vertical(0, g1, g2, g3, g4)

    def horizontal_recursion(self, r, g1, g2, g3, g4, g5):
        r_1 = g2.coordinates
        r_2 = g3.coordinates
        out1 = self.hgp_begin_horizontal(g1, g3, g4, g5)
        out2 = (r_1[r] - r_2[r]) * self.hgp_begin_horizontal(g2, g3, g4, g5)
        return out1 + out2

    def hgp_begin_vertical(self, m, g1, g2, g3, g4):
        l_1 = g1.integral_exponents
        l_3 = g3.integral_exponents
        if l_1[0] > 0:
            recursive_array = self.hgp_vertical_factory(0, g1, g2, g3, g4)
            return self.vertical_recursion(0, m, *recursive_array)
        if l_3[0] > 0:
            recursive_array = self.hgp_vertical_factory(0, g3, g4, g1, g2)
            return self.vertical_recursion(0, m, *recursive_array)
        if l_1[1] > 0:
            recursive_array = self.hgp_vertical_factory(1, g1, g2, g3, g4)
            return self.vertical_recursion(1, m, *recursive_array)
        if l_3[1] > 0:
            recursive_array = self.hgp_vertical_factory(1, g3, g4, g1, g2)
            return self.vertical_recursion(1, m, *recursive_array)
        if l_1[2] > 0:
            recursive_array = self.hgp_vertical_factory(2, g1, g2, g3, g4)
            return self.vertical_recursion(2, m, *recursive_array)
        if l_3[2] > 0:
            recursive_array = self.hgp_vertical_factory(2, g3, g4, g1, g2)
            return self.vertical_recursion(2, m, *recursive_array)
        else:
            return self.hgp_end(m, g1, g2, g3, g4)

    def vertical_recursion(self, r, m, g1, g2, g3, g4, g5, g6):
        a_1 = g1.exponent
        a_2 = g2.exponent
        a_3 = g3.exponent
        a_4 = g4.exponent
        a_5 = a_1 + a_2
        a_6 = a_3 + a_4
        rho = (a_5 * a_6) / (a_5 + a_6)

        r_1 = g1.coordinates
        r_2 = g2.coordinates
        r_3 = g3.coordinates
        r_4 = g4.coordinates
        r_5 = Vector.gaussian_product(a_1, r_1, a_2, r_2)
        r_6 = Vector.gaussian_product(a_3, r_3, a_4, r_4)
        r_7 = Vector.gaussian_product(a_5, r_5, a_6, r_6)

        out1 = (r_5[r] - r_1[r]) * self.hgp_begin_vertical(m, g1, g2, g3, g4)
        out2 = (r_7[r] - r_5[r]) * self.hgp_begin_vertical((m+1), g1, g2, g3, g4)
        if g5.integral_exponents[r] >= 0:
            out3 = self.os_count(g1.integral_exponents[r]) * (1 / (2 * a_5)) * self.hgp_begin_vertical(m, g5, g2, g3, g4)
            out4 = - self.os_count(g1.integral_exponents[r]) * (rho / (2 * a_5**2)) * self.hgp_begin_vertical((m+1), g5, g2, g3, g4)
        else:
            out3 = 0
            out4 = 0
        if g6.integral_exponents[r] >= 0:
            out5 = self.os_count(g3.integral_exponents[r]) * (1 / (2*(a_5 + a_6))) * self.hgp_begin_vertical((m+1), g1, g2, g6, g4)
        else:
            out5 = 0
        return out1 + out2 + out3 + out4 + out5

    @staticmethod
    def os_count(i):
        if i == 0:
            return 1
        else:
            return i

    @staticmethod
    def hgp_horizontal_factory(r, g1, g2, g3, g4):
        d_1 = g1.contraction
        d_2 = g2.contraction
        a_1 = g1.exponent
        a_2 = g2.exponent
        r_1 = g1.coordinates
        r_2 = g2.coordinates
        l_1 = g1.integral_exponents
        l_2 = g2.integral_exponents
        if r == 0:
            g1xa1 = PrimitiveBasis(d_1, a_1, r_1, (l_1[0] + 1, l_1[1], l_1[2]))
            g2xm1 = PrimitiveBasis(d_2, a_2, r_2, (l_2[0] - 1, l_2[1], l_2[2]))
            return g1xa1, g1, g2xm1, g3, g4
        elif r == 1:
            g1ya1 = PrimitiveBasis(d_1, a_1, r_1, (l_1[0], l_1[1] + 1, l_1[2]))
            g2ym1 = PrimitiveBasis(d_2, a_2, r_2, (l_2[0], l_2[1] - 1, l_2[2]))
            return g1ya1, g1, g2ym1, g3, g4
        elif r == 2:
            g1za1 = PrimitiveBasis(d_1, a_1, r_1, (l_1[0], l_1[1], l_1[2] + 1))
            g2zm1 = PrimitiveBasis(d_2, a_2, r_2, (l_2[0], l_2[1], l_2[2] - 1))
            return g1za1, g1, g2zm1, g3, g4

    @staticmethod
    def hgp_vertical_factory(r, g1, g2, g3, g4):
        d_1 = g1.contraction
        d_3 = g3.contraction
        a_1 = g1.exponent
        a_3 = g3.exponent
        r_1 = g1.coordinates
        r_3 = g3.coordinates
        l_1 = g1.integral_exponents
        l_3 = g3.integral_exponents
        if r == 0:
            g1xm1 = PrimitiveBasis(d_1, a_1, r_1, (l_1[0] - 1, l_1[1], l_1[2]))
            g1xm2 = PrimitiveBasis(d_1, a_1, r_1, (l_1[0] - 2, l_1[1], l_1[2]))
            g3xm1 = PrimitiveBasis(d_3, a_3, r_3, (l_3[0] - 1, l_3[1], l_3[2]))
            return g1xm1, g2, g3, g4, g1xm2, g3xm1
        elif r == 1:
            g1ym1 = PrimitiveBasis(d_1, a_1, r_1, (l_1[0], l_1[1] - 1, l_1[2]))
            g1ym2 = PrimitiveBasis(d_1, a_1, r_1, (l_1[0], l_1[1] - 2, l_1[2]))
            g3ym1 = PrimitiveBasis(d_3, a_3, r_3, (l_3[0], l_3[1] - 1, l_3[2]))
            return g1ym1, g2, g3, g4, g1ym2, g3ym1
        elif r == 2:
            g1zm1 = PrimitiveBasis(d_1, a_1, r_1, (l_1[0], l_1[1], l_1[2] - 1))
            g1zm2 = PrimitiveBasis(d_1, a_1, r_1, (l_1[0], l_1[1], l_1[2] - 2))
            g3zm1 = PrimitiveBasis(d_3, a_3, r_3, (l_3[0], l_3[1], l_3[2] - 1))
            return g1zm1, g2, g3, g4, g1zm2, g3zm1

    def hgp_end(self, m, g1, g2, g3, g4):
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
            r_5 = Vector.gaussian_product(a_1, r_1, a_2, r_2)
            r_6 = Vector.gaussian_product(a_3, r_3, a_4, r_4)
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
