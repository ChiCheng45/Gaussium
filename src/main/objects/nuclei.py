import numpy as np

class Nuclei:

    def __init__(self, array):
        self.name = array[0]
        self.charge = float(array[1])
        self.mass = float(array[2])
        self.coordinates = np.matrix([[float(array[3])], [float(array[4])], [float(array[5])]])

    def get_name(self):
        return self.name

    def get_charge(self):
        return self.charge

    def get_mass(self):
        return self.mass

    def get_coordinates(self):
        return self.coordinates
