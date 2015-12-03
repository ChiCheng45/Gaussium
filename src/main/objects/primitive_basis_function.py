from math import pi
from math import factorial as fac


class PrimitiveBasis:
    def __init__(self, contraction, exponent, coordinates, integral_exponents):
        self.contraction = contraction
        self.exponent = exponent
        self.coordinates = coordinates
        self.integral_exponents = integral_exponents

    def normalisation(self):
        normalisation = (((2 * self.exponent) / pi) ** (3 / 4)) * (((((8 * self.exponent) ** (
            self.integral_exponents[0] + self.integral_exponents[1] + self.integral_exponents[2])) * fac(
            self.integral_exponents[0]) * fac(self.integral_exponents[1]) * fac(self.integral_exponents[2])) / (fac(
            2 * self.integral_exponents[0]) * fac(2 * self.integral_exponents[1]) * fac(
            2 * self.integral_exponents[2]))) ** (1 / 2))
        return normalisation
