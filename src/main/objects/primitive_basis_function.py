class PrimitiveBasis:

    def __init__(self, orbital_type, contraction, exponent):
        self.orbital_type = orbital_type
        self.contraction = contraction
        self.exponent = exponent

    def get_orbital_type(self):
        return self.orbital_type

    def get_contraction(self):
        return self.contraction

    def get_exponent(self):
        return self.exponent
