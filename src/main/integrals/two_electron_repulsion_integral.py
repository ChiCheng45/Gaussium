from src.main.common import Binomial, Vector, BoysFunction
from math import factorial as fac
from math import sqrt, pi, exp
from src.main.objects import PrimitiveBasis


class ElectronRepulsionIntegral:

    @staticmethod
    def gaussian_product_factor(a_1, a_2, a_3, a_4, a_p, a_q, r_a, r_b, r_c, r_d):
        r_ab = Vector.distance(r_a, r_b)
        r_cd = Vector.distance(r_c, r_d)
        ans = ((2 * pi**2) / (a_p * a_q)) * sqrt(pi / (a_p + a_q)) * exp(- ((a_1 * a_2 * r_ab**2) / a_p) - ((a_3 * a_4 * r_cd**2) / a_q))
        return ans

    @staticmethod
    def sigma_function(l, l_1, l_2, a, b, r, g):
        ans = Binomial.calculate_coefficient(l, l_1, l_2, a, b) * ((fac(l) * g**(r - l)) / (fac(r) * fac(l - 2*r)))
        return ans

    @staticmethod
    def b_function(l, ll, r, rr, i, l_1, l_2, a_x, b_x, p_x, g_1, l_3, l_4, c_x, d_x, q_x, g_2):
        pa_x = p_x - a_x
        pb_x = p_x - b_x
        qc_x = q_x - c_x
        qd_x = q_x - d_x
        c_x = p_x - q_x
        delta = (1/(4*g_1)) + (1/(4*g_2))
        out1 = (-1)**l * ElectronRepulsionIntegral.sigma_function(l, l_1, l_2, pa_x, pb_x, r, g_1) * ElectronRepulsionIntegral.sigma_function(ll, l_3, l_4, qc_x, qd_x, rr, g_2)
        out2_num = (-1)**i * (2 * delta)**(2 * (r + rr)) * fac(l + ll - 2*r - 2*rr) * delta**i * c_x**(l + ll - 2*(r + rr + i))
        out2_den = (4 * delta)**(l + ll) * fac(i) * fac(l + ll - 2*(r + rr + i))
        ans = out1 * (out2_num / out2_den)
        return ans

    @staticmethod
    def integral(g1, g2, g3, g4):
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
        r_56 = Vector.distance(r_5, r_6)

        l_1 = g1.integral_exponents
        l_2 = g2.integral_exponents
        l_3 = g3.integral_exponents
        l_4 = g4.integral_exponents

        delta = (1/(4*a_5)) + (1/(4*a_6))

        out1 = 0
        for l in range(l_1[0] + l_2[0] + 1):
            for r in range(int(l/2) + 1):
                for ll in range(l_3[0] + l_4[0] + 1):
                    for rr in range(int(ll/2) + 1):
                        for i in range(int((l + ll - 2*r - 2*rr) / 2) + 1):
                            out2 = ElectronRepulsionIntegral.b_function(l, ll, r, rr, i, l_1[0], l_2[0], r_1[0], r_2[0], r_5[0], a_5, l_3[0], l_4[0], r_3[0], r_4[0], r_6[0], a_6)
                            for m in range(l_1[1] + l_2[1] + 1):
                                for s in range(int(m / 2) + 1):
                                    for mm in range(l_3[1] + l_4[1] + 1):
                                        for ss in range(int(mm/2) + 1):
                                            for j in range(int((m + mm - 2*s - 2*ss) / 2) + 1):
                                                out3 = ElectronRepulsionIntegral.b_function(m, mm, s, ss, j, l_1[1], l_2[1], r_1[1], r_2[1], r_5[1], a_5, l_3[1], l_4[1], r_3[1], r_4[1], r_6[1], a_6)
                                                for n in range(l_1[2] + l_2[2] + 1):
                                                    for t in range(int(n/2) + 1):
                                                        for nn in range(l_3[2] + l_4[2] + 1):
                                                            for tt in range(int(nn/2) + 1):
                                                                for k in range(int((n + nn - 2*t - 2*tt) / 2) + 1):
                                                                    out4 = ElectronRepulsionIntegral.b_function(n, nn, t, tt, k, l_1[2], l_2[2], r_1[2], r_2[2], r_5[2], a_5, l_3[2], l_4[2], r_3[2], r_4[2], r_6[2], a_6)
                                                                    v = l + ll + m + mm + n + nn - 2*(r + rr + s + ss + t + tt) - (i + j + k)
                                                                    out1 += out2 * out3 * out4 * BoysFunction.function(v, (r_56**2 / (4 * delta)))
        out1 *= ElectronRepulsionIntegral.gaussian_product_factor(a_1, a_2, a_3, a_4, a_5, a_6, r_1, r_2, r_3, r_4)
        return out1


class ObaraSaika:

    @staticmethod
    def os_gaussian_factory(g1, g2, g3, g4, coord):
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

        if coord == 0:
            g1x1 = PrimitiveBasis(d_1, a_1, r_1, (l_1[0] - 1, l_1[1], l_1[2]))
            g2x1 = PrimitiveBasis(d_2, a_2, r_2, (l_2[0] - 1, l_2[1], l_2[2]))
            g3x1 = PrimitiveBasis(d_3, a_3, r_3, (l_3[0] - 1, l_3[1], l_3[2]))
            g4x1 = PrimitiveBasis(d_4, a_4, r_4, (l_4[0] - 1, l_4[1], l_4[2]))
            return tuple((g1x1, g2, g3, g4, g2x1, g3x1, g4x1))
        elif coord == 1:
            g1y1 = PrimitiveBasis(d_1, a_1, r_1, (l_1[0], l_1[1] - 1, l_1[2]))
            g2y1 = PrimitiveBasis(d_2, a_2, r_2, (l_2[0], l_2[1] - 1, l_2[2]))
            g3y1 = PrimitiveBasis(d_3, a_3, r_3, (l_3[0], l_3[1] - 1, l_3[2]))
            g4y1 = PrimitiveBasis(d_4, a_4, r_4, (l_4[0], l_4[1] - 1, l_4[2]))
            return tuple((g1y1, g2, g3, g4, g2y1, g3y1, g4y1))
        elif coord == 2:
            g1z1 = PrimitiveBasis(d_1, a_1, r_1, (l_1[0], l_1[1], l_1[2] - 1))
            g2z1 = PrimitiveBasis(d_2, a_2, r_2, (l_2[0], l_2[1], l_2[2] - 1))
            g3z1 = PrimitiveBasis(d_3, a_3, r_3, (l_3[0], l_3[1], l_3[2] - 1))
            g4z1 = PrimitiveBasis(d_4, a_4, r_4, (l_4[0], l_4[1], l_4[2] - 1))
            return tuple((g1z1, g2, g3, g4, g2z1, g3z1, g4z1))

    @staticmethod
    def os_set(g1, g2, g3, g4):
        return ObaraSaika.os_begin(g1, g2, g3, g4, 0)

    @staticmethod
    def os_begin(g1, g2, g3, g4, m):
        l_1 = g1.integral_exponents
        l_2 = g2.integral_exponents
        l_3 = g3.integral_exponents
        l_4 = g4.integral_exponents

        if l_1[0] + l_1[1] + l_1[2] + l_2[0] + l_2[1] + l_2[2] + l_3[0] + l_3[1] + l_3[2] + l_4[0] + l_4[1] + l_4[2] == 0:
            return ObaraSaika.os_end(g1, g2, g3, g4, m)
        else:
            if l_1[0] > 0:
                recursive_array = ObaraSaika.os_gaussian_factory(g1, g2, g3, g4, 0)
                return ObaraSaika.os_recursive(recursive_array, 0, m)
            if l_1[1] > 0:
                recursive_array = ObaraSaika.os_gaussian_factory(g1, g2, g3, g4, 1)
                return ObaraSaika.os_recursive(recursive_array, 1, m)
            if l_1[2] > 0:
                recursive_array = ObaraSaika.os_gaussian_factory(g1, g2, g3, g4, 2)
                return ObaraSaika.os_recursive(recursive_array, 2, m)
            if l_2[0] > 0:
                recursive_array = ObaraSaika.os_gaussian_factory(g2, g1, g3, g4, 0)
                return ObaraSaika.os_recursive(recursive_array, 0, m)
            if l_2[1] > 0:
                recursive_array = ObaraSaika.os_gaussian_factory(g2, g1, g3, g4, 1)
                return ObaraSaika.os_recursive(recursive_array, 1, m)
            if l_2[2] > 0:
                recursive_array = ObaraSaika.os_gaussian_factory(g2, g1, g3, g4, 2)
                return ObaraSaika.os_recursive(recursive_array, 2, m)
            if l_3[0] > 0:
                recursive_array = ObaraSaika.os_gaussian_factory(g3, g4, g1, g2, 0)
                return ObaraSaika.os_recursive(recursive_array, 0, m)
            if l_3[1] > 0:
                recursive_array = ObaraSaika.os_gaussian_factory(g3, g4, g1, g2, 1)
                return ObaraSaika.os_recursive(recursive_array, 1, m)
            if l_3[2] > 0:
                recursive_array = ObaraSaika.os_gaussian_factory(g3, g4, g1, g2, 2)
                return ObaraSaika.os_recursive(recursive_array, 2, m)
            if l_4[0] > 0:
                recursive_array = ObaraSaika.os_gaussian_factory(g4, g3, g1, g2, 0)
                return ObaraSaika.os_recursive(recursive_array, 0, m)
            if l_4[1] > 0:
                recursive_array = ObaraSaika.os_gaussian_factory(g4, g3, g1, g2, 1)
                return ObaraSaika.os_recursive(recursive_array, 1, m)
            if l_4[2] > 0:
                recursive_array = ObaraSaika.os_gaussian_factory(g4, g3, g1, g2, 2)
                return ObaraSaika.os_recursive(recursive_array, 2, m)

    @staticmethod
    def os_recursive(g, coord, m):
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

        l_2 = g[4].integral_exponents[coord]
        l_3 = g[5].integral_exponents[coord]
        l_4 = g[6].integral_exponents[coord]

        out1 = (r_5[coord] - r_1[coord]) * ObaraSaika.os_begin(g[0], g[1], g[2], g[3], m)
        out2 = (r_7[coord] - r_5[coord]) * ObaraSaika.os_begin(g[0], g[1], g[2], g[3], (m+1))
        if l_2 >= 0:
            out5 = (1 / (2 * a_5)) * ObaraSaika.os_begin(g[0], g[4], g[2], g[3], m)
            out6 = - (rho / (2 * a_5**2)) * ObaraSaika.os_begin(g[0], g[4], g[2], g[3], (m+1))
        else:
            out5 = 0
            out6 = 0
        if l_3 >= 0:
            out7 = (1 / (2*(a_5 + a_6))) * ObaraSaika.os_begin(g[0], g[1], g[5], g[3], (m+1))
        else:
            out7 = 0
        if l_4 >= 0:
            out8 = (1 / (2*(a_5 + a_6))) * ObaraSaika.os_begin(g[0], g[1], g[2], g[6], (m+1))
        else:
            out8 = 0
        return out1 + out2 + out5 + out6 + out7 + out8

    @staticmethod
    def os_end(g1, g2, g3, g4, m):
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
        return out1 * out2 * out3
