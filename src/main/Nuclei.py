class Nuclei:

    def __init__(self, array):
        self.name = array[0]
        self.charge = float(array[1])
        self.mass = float(array[2])
        self.x = float(array[3])
        self.y = float(array[4])
        self.z = float(array[5])

    def get_name(self):
        return self.name

    def get_charge(self):
        return self.charge

    def get_mass(self):
        return self.mass

    def get_x(self):
        return self.x

    def get_y(self):
        return self.y

    def get_z(self):
        return self.z
