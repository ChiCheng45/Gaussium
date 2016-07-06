from src.main.common import create_quaternion
from src.main.common import quaternion_rotation
from math import pi


class RotationSymmetry:

    def __init__(self, fold, vector):
        self.fold = fold
        self.vector = vector

    def operate(self, coordinate):
        angle = 2 * pi / self.fold
        quaternion = create_quaternion(self.vector, angle)
        return quaternion_rotation(quaternion, coordinate)

    def int_operate(self, coordinate):
        angle = 2 * pi / self.fold
        quaternion = create_quaternion(self.vector, angle)
        x, y, z = quaternion_rotation(quaternion, coordinate)
        return int(round(x, 1)), int(round(y, 1)), int(round(z, 1))

    @property
    def symmetry_operation(self):
        if -1e-3 <= self.vector[0] <= 1e-3:
            i = 0.0
        else:
            i = round(self.vector[0], 3)
        if -1e-3 <= self.vector[1] <= 1e-3:
            j = 0.0
        else:
            j = round(self.vector[1], 3)
        if -1e-3 <= self.vector[2] <= 1e-3:
            k = 0.0
        else:
            k = round(self.vector[2], 3)
        return 'C_{' + str(round(self.fold, 3)) + '}' + str((i, j, k))
