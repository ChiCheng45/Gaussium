class Basis:

    def __init__(self, name, coordinates, primitive_gaussian_array):
        self.name = name
        self.coordinates = coordinates
        self.primitive_gaussian_array = primitive_gaussian_array

    def get_name(self):
        return self.name

    def get_coordinates(self):
        return self.coordinates

    def get_primitive_gaussian_array(self):
        return self.primitive_gaussian_array
