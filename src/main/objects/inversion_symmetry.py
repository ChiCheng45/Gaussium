class InversionSymmetry:

    def __init__(self):
        self.symmetry_operation = 'i'

    def operate(self, coordinate):
        return -coordinate[0], -coordinate[1], -coordinate[2]
