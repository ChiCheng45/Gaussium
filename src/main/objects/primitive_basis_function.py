from math import pi
from math import factorial as fac


class PrimitiveBasis:
    def __init__(self, contraction, exponent, coordinates, integral_exponents):
        self.contraction = contraction
        self.exponent = exponent
        self.coordinates = coordinates
        self.integral_exponents = integral_exponents
        self.normalisation = (((2 * exponent) / pi) ** (3 / 4)) * (((((8 * exponent) ** (
            integral_exponents[0] + integral_exponents[1] + integral_exponents[2])) * fac(integral_exponents[0]) * fac(
            integral_exponents[1]) * fac(integral_exponents[2])) / (fac(2 * integral_exponents[0]) * fac(
            2 * integral_exponents[1]) * fac(2 * integral_exponents[2]))) ** (1 / 2))
