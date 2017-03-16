from src.common import create_householder_matrix
from src.common import householder_matrix_reflection


class ReflectionSymmetry:
    """Reflection symmetry operator.

    This class can only create reflection planes that go through the origin.

    Attributes
    ----------
    vector : Tuple[float, float, float]
        A vector that is orthogonal to the reflection plane.

    """
    def __init__(self, vector):
        self.vector = vector

    def operate(self, coordinate):
        """Reflects a point through the reflection plane.

        Parameters
        ----------
        coordinate : Tuple[float, float, float]

        Returns
        -------
        : Tuple[float, float, float]

        """
        householder_matrix = create_householder_matrix(self.vector)
        return householder_matrix_reflection(coordinate, householder_matrix)

    def int_operate(self, coordinate):
        """Returns the reflection symmetry operation on a coordinate but returns a tuple of ints.

        This method is useful for seeing how the symmetry operation affects the sign of a gaussian functions integral
        exponents.

        Parameters
        ----------
        coordinate : Tuple[float, float, float]

        Returns
        -------
        : Tuple[int, int, int]

        """
        householder_matrix = create_householder_matrix(self.vector)
        x, y, z = householder_matrix_reflection(coordinate, householder_matrix)
        return int(round(x, 1)), int(round(y, 1)), int(round(z, 1))

    @property
    def symmetry_operation(self):
        """Returns the reflection symmetry label.

        Returns
        -------
        : str

        """
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
