import numpy as np


class Basis:

    def __init__(self, name, coordinates, orbital_type, array_of_coefficients):
        self.name = name
        self.coordinates = coordinates
        self.orbital_type = orbital_type
        self.array_of_coefficients = array_of_coefficients

    def get_name(self):
        return self.name

    def get_coordinates(self):
        return self.coordinates

    def get_orbital_type(self):
        return self.orbital_type

    def get_array_of_coefficients(self):
        return self.array_of_coefficients
