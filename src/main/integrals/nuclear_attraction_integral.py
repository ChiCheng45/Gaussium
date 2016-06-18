from src.main.integrals import Binomial
from src.main.integrals import BoysFunction
from src.main.common import gaussian_product_coordinate
from src.main.common import coordinate_distance
from src.main.common import vector_minus
from math import factorial as fac
from math import exp, pi


class NuclearAttraction:

    @staticmethod
    def a_function(l, r, i, l_1, l_2, pa, pb, pc, g):
        e = 1 / (4*g)
        out1 = Binomial.coefficient(l, l_1, l_2, pa, pb)
        out2 = (-1)**i * fac(l) * pc**(l - 2*r - 2*i) * e**(r + i)
        out3 = fac(r) * fac(i) * fac(l - 2*r - 2*i)
        ans = (-1)**l * out1 * (out2/out3)
        return ans

    @classmethod
    def integral(cls, gaussian_1, gaussian_2, nuclei):
        a_1 = gaussian_1.exponent
        a_2 = gaussian_2.exponent
        l_1 = gaussian_1.integral_exponents
        l_2 = gaussian_2.integral_exponents

        r_a = gaussian_1.coordinates
        r_b = gaussian_2.coordinates
        r_c = nuclei.coordinates
        r_p = gaussian_product_coordinate(a_1, r_a, a_2, r_b)

        r_ab = coordinate_distance(r_a, r_b)
        r_pc = coordinate_distance(r_p, r_c)

        r_p_a = vector_minus(r_p, r_a)
        r_p_b = vector_minus(r_p, r_b)
        r_p_c = vector_minus(r_p, r_c)

        g = a_1 + a_2

        ans = 0
        for l in range(l_1[0] + l_2[0] + 1):
            for r in range(int(l/2) + 1):
                for i in range(int((l - 2*r) / 2) + 1):
                    out1 = cls.a_function(l, r, i, l_1[0], l_2[0], r_p_a[0], r_p_b[0], r_p_c[0], g)
                    for m in range(l_1[1] + l_2[1] + 1):
                        for s in range(int(m/2) + 1):
                            for j in range(int((m - 2*s) / 2) + 1):
                                out2 = cls.a_function(m, s, j, l_1[1], l_2[1], r_p_a[1], r_p_b[1], r_p_c[1], g)
                                for n in range(l_1[2] + l_2[2] + 1):
                                    for t in range(int(n/2) + 1):
                                        for k in range(int((n - 2*t) / 2) + 1):
                                            out3 = cls.a_function(n, t, k, l_1[2], l_2[2], r_p_a[2], r_p_b[2], r_p_c[2], g)
                                            v = (l + m + n) - 2*(r + s + t) - (i + j + k)
                                            out4 = BoysFunction.calculate(v, g * r_pc**2)
                                            out5 = out1 * out2 * out3 * out4
                                            ans += out5
        ans *= ((2 * pi) / g) * exp(- (a_1 * a_2 * r_ab**2) / g)
        return ans
