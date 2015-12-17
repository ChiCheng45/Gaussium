from math import pi, sqrt
from math import factorial as fac


class PrimitiveBasis:
    def __init__(self, contraction, exponent, coordinates, integral_exponents):
        self.contraction = contraction
        self.exponent = exponent
        self.coordinates = coordinates
        self.integral_exponents = integral_exponents

    def normalisation(self):
        l = self.integral_exponents
        out1 = ((2 * self.exponent) / pi) ** (3 / 4)
        out2 = (8 * self.exponent) ** (l.l + l.m + l.n) * fac(l.l) * fac(l.m) * fac(l.n)
        out3 = fac(2 * l.l) * fac(2 * l.m) * fac(2 * l.n)
        normalisation = out1 * sqrt(out2 / out3)
        return normalisation
