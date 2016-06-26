from src.main.common import create_householder_matrix
from src.main.common import householder_matrix_reflection


class ReflectionSymmetry:

    def __init__(self, vector):
        self.vector = vector

    def operate(self, coordinate):
        householder_matrix = create_householder_matrix(self.vector)
        return householder_matrix_reflection(coordinate, householder_matrix)

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
        return 'sigma' + str((i, j, k))
