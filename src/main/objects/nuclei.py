class Nuclei:

    def __init__(self, array):
        self.element = array[0]
        self.charge = float(array[1])
        self.mass = float(array[2])
        self.coordinates = float(array[3]), float(array[4]), float(array[5])
