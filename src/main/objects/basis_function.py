class Basis:

    def __init__(self, name, x, y, z, orbital_type, array_of_coefficients):
        self.name = name
        self.x = x
        self.y = y
        self.z = z
        self.orbital_type = orbital_type
        self.array_of_coefficients = array_of_coefficients

    def get_name(self):
        return self.name

    def get_x(self):
        return self.x

    def get_y(self):
        return self.y

    def get_z(self):
        return self.z

    def get_orbital_type(self):
        return self.orbital_type

    def get_array_of_coefficients(self):
        return self.array_of_coefficients
