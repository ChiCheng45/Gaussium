import itertools
from math import factorial as fac
from math import sqrt, pi, exp
from gaussium.common import coordinate_distance
from gaussium.factory import gaussian_product
from gaussium.integrals import binomial_coefficient
from gaussium.integrals import boys_function


class ElectronRepulsion:

    def __init__(self):
        self.end_dict = {}

    def integrate(self, basis_i, basis_j, basis_k, basis_l):
        primitives_i = basis_i.primitive_gaussian_array
        primitives_j = basis_j.primitive_gaussian_array
        primitives_k = basis_k.primitive_gaussian_array
        primitives_l = basis_l.primitive_gaussian_array

        l_1 = basis_i.integral_exponents
        l_2 = basis_j.integral_exponents
        l_3 = basis_k.integral_exponents
        l_4 = basis_l.integral_exponents

        r_1 = basis_i.coordinates
        r_2 = basis_j.coordinates
        r_3 = basis_k.coordinates
        r_4 = basis_l.coordinates

        n_i = basis_i.normalisation
        n_j = basis_j.normalisation
        n_k = basis_k.normalisation
        n_l = basis_l.normalisation
        norm = n_i * n_j * n_k * n_l

        ans = 0.0
        for g1, g2, g3, g4 in itertools.product(primitives_i, primitives_j, primitives_k, primitives_l):
            self.end_dict = {}
            c_1 = g1.contraction
            c_2 = g2.contraction
            c_3 = g3.contraction
            c_4 = g4.contraction
            n_1 = g1.normalisation
            n_2 = g2.normalisation
            n_3 = g3.normalisation
            n_4 = g4.normalisation
            contraction = c_1 * c_2 * c_3 * c_4 * n_1 * n_2 * n_3 * n_4 * norm

            g5 = gaussian_product(g1, g2)
            g6 = gaussian_product(g3, g4)

            a_1 = g1.exponent
            a_2 = g2.exponent
            a_3 = g3.exponent
            a_4 = g4.exponent
            a_5 = g5.exponent
            a_6 = g6.exponent

            r_5 = g5.coordinates
            r_6 = g6.coordinates
            r_12 = coordinate_distance(r_1, r_2)
            r_34 = coordinate_distance(r_3, r_4)
            r_56 = coordinate_distance(r_5, r_6)

            l_5 = g5.integral_exponents
            l_6 = g6.integral_exponents

            delta = (1/(4*a_5)) + (1/(4*a_6))

            out_5 = 0
            for l in range(l_5[0] + 1):
                for r in range(int(l/2) + 1):
                    for ll in range(l_6[0] + 1):
                        for rr in range(int(ll/2) + 1):
                            for i in range(int((l + ll - 2*r - 2*rr) / 2) + 1):
                                out1 = self.b_function(l, ll, r, rr, i, l_1[0], l_2[0], r_1[0], r_2[0], r_5[0], a_5, l_3[0], l_4[0], r_3[0], r_4[0], r_6[0], a_6)
                                for m in range(l_5[1] + 1):
                                    for s in range(int(m/2) + 1):
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
                                                                        out_5 += out1 * out2 * out3 * out4
            out_5 *= self.gaussian_product_factor(a_1, a_2, a_3, a_4, a_5, a_6, r_12, r_34)
            ans += contraction * out_5
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
        return binomial_coefficient(l, l_1, l_2, a, b) * ((fac(l) * g**(r - l)) / (fac(r) * fac(l - 2*r)))

    def gaussian_product_factor(self, a_1, a_2, a_3, a_4, a_p, a_q, r_ab, r_cd):
        ans = ((2 * pi**2) / (a_p * a_q)) * sqrt(pi / (a_p + a_q)) * exp(- ((a_1 * a_2 * r_ab**2) / a_p) - ((a_3 * a_4 * r_cd**2) / a_q))
        return ans
