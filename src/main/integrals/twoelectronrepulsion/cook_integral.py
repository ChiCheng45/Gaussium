from src.main.common import Binomial, Vector, BoysFunction
from src.main.objects import PrimitiveBasisFactory
from math import factorial as fac
from math import sqrt, pi, exp


class ElectronRepulsionIntegral:

    def __init__(self):
        self.end_dict = {}

    def integral(self, g1, g2, g3, g4):
        self.end_dict = {}
        boys_function = BoysFunction.function

        g5 = PrimitiveBasisFactory.gaussian_product(g1, g2)
        g6 = PrimitiveBasisFactory.gaussian_product(g3, g4)

        a_1 = g1.exponent
        a_2 = g2.exponent
        a_3 = g3.exponent
        a_4 = g4.exponent
        a_5 = g5.exponent
        a_6 = g6.exponent

        r_1 = g1.coordinates
        r_2 = g2.coordinates
        r_3 = g3.coordinates
        r_4 = g4.coordinates
        r_5 = g5.coordinates
        r_6 = g6.coordinates
        r_12 = Vector.distance(r_1, r_2)
        r_34 = Vector.distance(r_3, r_4)
        r_56 = Vector.distance(r_5, r_6)

        l_1 = g1.integral_exponents
        l_2 = g2.integral_exponents
        l_3 = g3.integral_exponents
        l_4 = g4.integral_exponents
        l_5 = g5.integral_exponents
        l_6 = g6.integral_exponents

        delta = (1/(4*a_5)) + (1/(4*a_6))

        ans = 0
        for l in range(l_5[0] + 1):
            for r in range(int(l/2) + 1):
                for ll in range(l_6[0] + 1):
                    for rr in range(int(ll/2) + 1):
                        for i in range(int((l + ll - 2*r - 2*rr) / 2) + 1):
                            out1 = self.b_function(l, ll, r, rr, i, l_1[0], l_2[0], r_1[0], r_2[0], r_5[0], a_5, l_3[0], l_4[0], r_3[0], r_4[0], r_6[0], a_6)
                            for m in range(l_5[1] + 1):
                                for s in range(int(m / 2) + 1):
                                    for mm in range(l_6[1] + 1):
                                        for ss in range(int(mm/2) + 1):
                                            for j in range(int((m + mm - 2*s - 2*ss) / 2) + 1):
                                                out2 = self.b_function(m, mm, s, ss, j, l_1[1], l_2[1], r_1[1], r_2[1], r_5[1], a_5, l_3[1], l_4[1], r_3[1], r_4[1], r_6[1], a_6)
                                                for n in range(l_5[2] + 1):
                                                    for t in range(int(n/2) + 1):
                                                        for nn in range(l_6[2] + 1):
                                                            for tt in range(int(nn/2) + 1):
                                                                for k in range(int((n + nn - 2*t - 2*tt) / 2) + 1):
                                                                    out3 = self.b_function(n, nn, t, tt, k, l_1[2], l_2[2], r_1[2], r_2[2], r_5[2], a_5, l_3[2], l_4[2], r_3[2], r_4[2], r_6[2], a_6)
                                                                    v = l + ll + m + mm + n + nn - 2*(r + rr + s + ss + t + tt) - (i + j + k)
                                                                    if v in self.end_dict:
                                                                        out4 = self.end_dict[v]
                                                                    else:
                                                                        out4 = boys_function(v, (r_56**2 / (4 * delta)))
                                                                        self.end_dict[v] = out4
                                                                    ans += out1 * out2 * out3 * out4
        ans *= self.gaussian_product_factor(a_1, a_2, a_3, a_4, a_5, a_6, r_12, r_34)
        return ans

    def b_function(self, l, ll, r, rr, i, l_1, l_2, a_x, b_x, p_x, g_1, l_3, l_4, c_x, d_x, q_x, g_2):
        sigma = self.sigma
        pa_x = p_x - a_x
        pb_x = p_x - b_x
        qc_x = q_x - c_x
        qd_x = q_x - d_x
        c_x = p_x - q_x
        delta = (1/(4*g_1)) + (1/(4*g_2))
        out1 = (-1)**l * sigma(l, l_1, l_2, pa_x, pb_x, r, g_1) * sigma(ll, l_3, l_4, qc_x, qd_x, rr, g_2)
        out2 = (-1)**i * (2 * delta)**(2 * (r + rr)) * fac(l + ll - 2*r - 2*rr) * delta**i * c_x**(l + ll - 2*(r + rr + i))
        out3 = (4 * delta)**(l + ll) * fac(i) * fac(l + ll - 2*(r + rr + i))
        ans = out1 * (out2 / out3)
        return ans

    def sigma(self, l, l_1, l_2, a, b, r, g):
        return Binomial.calculate_coefficient(l, l_1, l_2, a, b) * ((fac(l) * g**(r - l)) / (fac(r) * fac(l - 2*r)))

    def gaussian_product_factor(self, a_1, a_2, a_3, a_4, a_p, a_q, r_ab, r_cd):
        ans = ((2 * pi**2) / (a_p * a_q)) * sqrt(pi / (a_p + a_q)) * exp(- ((a_1 * a_2 * r_ab**2) / a_p) - ((a_3 * a_4 * r_cd**2) / a_q))
        return ans
