from math import pi, sqrt
from math import factorial as fac


class PrimitiveBasis:

    def __init__(self, contraction, exponent, coordinates, integral_exponents):
        self.contraction = contraction
        self.exponent = exponent
        self.coordinates = coordinates
        self.integral_exponents = integral_exponents
        self.normalisation_memo = None

    @property
    def normalisation(self):
        if self.normalisation_memo is None:
            x, y, z = self.integral_exponents
            out1 = ((2 * self.exponent) / pi) ** (3 / 4)
            out2 = (8 * self.exponent) ** (x + y + z) * fac(x) * fac(y) * fac(z)
            out3 = fac(2 * x) * fac(2 * y) * fac(2 * z)
            self.normalisation_memo = out1 * sqrt(out2 / out3)
        return self.normalisation_memo
