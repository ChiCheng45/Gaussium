from src.main.common import Binomial
from src.main.common import BoysFunction
from src.main.common import Vector
from math import factorial as fac
from math import exp, pi


class NuclearAttractionIntegral:

    @staticmethod
    def a_function(l, r, i, l_1, l_2, pa, pb, pc, g):
        e = 1 / (4*g)
        f_l = Binomial.calculate_coefficient(l, l_1, l_2, pa, pb)
        num = (-1)**i * fac(l) * pc**(l - 2*r - 2*i) * e**(r + i)
        dom = fac(r) * fac(i) * fac(l - 2*r - 2*i)
        out = (-1)**l * f_l * (num/dom)
        return out

    @classmethod
    def primitive_nuclear_attraction(cls, gaussian_1, gaussian_2, nuclei):
        a_1 = gaussian_1.exponent
        a_2 = gaussian_2.exponent
        l_1 = gaussian_1.integral_exponents
        l_2 = gaussian_2.integral_exponents

        r_a = gaussian_1.coordinates
        r_b = gaussian_2.coordinates
        r_c = nuclei.coordinates
        r_p = Vector.gaussian_product(a_1, r_a, a_2, r_b)

        r_ab = Vector.distance(r_a, r_b)
        r_pc = Vector.distance(r_p, r_c)

        r_p_a = Vector.minus(r_p, r_a)
        r_p_b = Vector.minus(r_p, r_b)
        r_p_c = Vector.minus(r_p, r_c)

        g = a_1 + a_2

        out1 = 0
        for l in range(l_1[0] + l_2[0] + 1):
            for r in range(int(l/2) + 1):
                for i in range(int((l - 2*r) / 2) + 1):
                    out2 = cls.a_function(l, r, i, l_1[0], l_2[0], r_p_a[0], r_p_b[0], r_p_c[0], g)
                    for m in range(l_1[1] + l_2[1] + 1):
                        for s in range(int(m/2) + 1):
                            for j in range(int((m - 2*s) / 2) + 1):
                                out3 = cls.a_function(m, s, j, l_1[1], l_2[1], r_p_a[1], r_p_b[1], r_p_c[1], g)
                                for n in range(l_1[2] + l_2[2] + 1):
                                    for t in range(int(n/2) + 1):
                                        for k in range(int((n - 2*t) / 2) + 1):
                                            out4 = cls.a_function(n, t, k, l_1[2], l_2[2], r_p_a[2], r_p_b[2], r_p_c[2], g)
                                            v = (l + m + n) - 2*(r + s + t) - (i + j + k)
                                            out5 = BoysFunction.function(v, g * r_pc**2)
                                            out6 = out2 * out3 * out4 * out5
                                            out1 += out6
        out1 *= ((2 * pi) / g) * exp(- (a_1 * a_2 * r_ab**2) / g)
        return out1
