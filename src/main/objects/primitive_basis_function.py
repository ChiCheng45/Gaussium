class PrimitiveBasis:

    def __init__(self, orbital_type, coefficients):
        self.orbital_type = orbital_type
        self.coefficients = coefficients

    def get_orbital_type(self):
        return self.orbital_type

    def get_coefficients(self):
        return self.coefficients
