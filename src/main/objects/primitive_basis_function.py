class PrimitiveBasis:

    def __init__(self, orbital_type, contraction, exponent, integral_exponents):
        self.orbital_type = orbital_type
        self.contraction = contraction
        self.exponent = exponent
        self.integral_exponents = integral_exponents

    def get_orbital_type(self):
        return self.orbital_type

    def get_contraction(self):
        return self.contraction

    def get_exponent(self):
        return self.exponent

    def get_integral_exponents(self):
        return self.integral_exponents
