class InversionSymmetry:

    def __init__(self):
        self.symmetry_operation = 'i'

    def operate(self, coordinate):
        return -coordinate[0], -coordinate[1], -coordinate[2]

    def int_operate(self, coordinate):
        x, y, z = coordinate
        return int(round(x, 1)), int(round(y, 1)), int(round(z, 1))
