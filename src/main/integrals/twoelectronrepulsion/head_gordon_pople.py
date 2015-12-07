from src.main.common import Vector, BoysFunction
from math import sqrt, pi, exp
from src.main.objects import PrimitiveBasis


class HeadGordonPople:

    def __init__(self):
        self.end_dict = {}

    def hgp_set(self, g1, g2, g3, g4):
        self.end_dict = {}
        return self.hgp_begin(g1, g2, g3, g4, 0)

    def hgp_begin(self, g1, g2, g3, g4, m):
        l_1 = g1.integral_exponents
        l_2 = g2.integral_exponents
        l_3 = g3.integral_exponents
        l_4 = g4.integral_exponents

        if l_1[0] + l_1[1] + l_1[2] + l_2[0] + l_2[1] + l_2[2] + l_3[0] + l_3[1] + l_3[2] + l_4[0] + l_4[1] + l_4[2] == 0:
            return self.hgp_end(g1, g2, g3, g4, m)

        if l_1[0] and l_2[0] > 0:
            if l_1[0] >= l_2[0]:
                recursive_array = self.hgp_horizontal_factory(g1, g2, g3, g4, 0)
                return self.horizontal_recursion(recursive_array, 0, m)
            else:
                recursive_array = self.hgp_horizontal_factory(g2, g1, g4, g3, 0)
                return self.horizontal_recursion(recursive_array, 0, m)
        if l_3[0] and l_4[0] > 0:
            if l_3[0] >= l_4[0]:
                recursive_array = self.hgp_horizontal_factory(g3, g4, g1, g2, 0)
                return self.horizontal_recursion(recursive_array, 0, m)
            else:
                recursive_array = self.hgp_horizontal_factory(g4, g3, g2, g1, 0)
                return self.horizontal_recursion(recursive_array, 0, m)
        if l_1[0] and l_3[0] == 0:
            recursive_array = self.hgp_vertical_factory(g1, g2, g3, g4, 0)
            return self.vertical_recursion(recursive_array, 0, m)
        elif l_1[0] and l_4[0] == 0:
            recursive_array = self.hgp_vertical_factory(g1, g2, g4, g3, 0)
            return self.vertical_recursion(recursive_array, 0, m)
        elif l_2[0] and l_3[0] == 0:
            recursive_array = self.hgp_vertical_factory(g2, g1, g3, g4, 0)
            return self.vertical_recursion(recursive_array, 0, m)
        elif l_2[0] and l_4[0] == 0:
            recursive_array = self.hgp_vertical_factory(g2, g1, g4, g3, 0)
            return self.vertical_recursion(recursive_array, 0, m)

    def horizontal_recursion(self, g, xyz, m):
        pass

    def vertical_recursion(self, g, xyz, m):
        pass

    @staticmethod
    def hgp_horizontal_factory(g1, g2, g3, g4, xyz):
        d_1 = g1.contraction
        d_2 = g2.contraction

        a_1 = g1.exponent
        a_2 = g2.exponent

        r_1 = g1.coordinates
        r_2 = g2.coordinates

        l_1 = g1.integral_exponents
        l_2 = g2.integral_exponents

        if xyz == 0:
            g1xa1 = PrimitiveBasis(d_1, a_1, r_1, (l_1[0] + 1, l_1[1], l_1[2]))
            g2x1 = PrimitiveBasis(d_2, a_2, r_2, (l_2[0] - 1, l_2[1], l_2[2]))
            return g1xa1, g1, g2x1, g3, g4
        elif xyz == 1:
            g1ya1 = PrimitiveBasis(d_1, a_1, r_1, (l_1[0], l_1[1] + 1, l_1[2]))
            g2y1 = PrimitiveBasis(d_2, a_2, r_2, (l_2[0], l_2[1] - 1, l_2[2]))
            return g1ya1, g1, g2y1, g3, g4
        elif xyz == 2:
            g1za1 = PrimitiveBasis(d_1, a_1, r_1, (l_1[0], l_1[1] + 1, l_1[2]))
            g2z1 = PrimitiveBasis(d_2, a_2, r_2, (l_2[0], l_2[1] - 1, l_2[2]))
            return g1za1, g1, g2z1, g3, g4

    @staticmethod
    def hgp_vertical_factory(g1, g2, g3, g4, xyz):
        return 0

    def hgp_end(self, g1, g2, g3, g4, m):
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
