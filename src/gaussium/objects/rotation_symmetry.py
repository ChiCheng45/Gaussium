from math import pi
from gaussium.common import create_quaternion
from gaussium.common import quaternion_rotation


class RotationSymmetry:
    """Rotation symmetry operator.

    This class can only create axis of rotation that go through
    the origin.

    Attributes
    ----------
    fold : float
        The fraction of 2 * pi angle a coordinate is rotated.
    vector : Tuple[float, float, float]
        The axis of rotation.

    """
    def __init__(self, fold, vector):
        self.fold = fold
        self.vector = vector

    def operate(self, coordinate):
        """Rotates a point around the axis of rotation for a angle of
        2 * pi / self.fold

        Parameters
        ----------
        coordinate : Tuple[float, float, float]

        Returns
        -------
        : Tuple[float, float, float]

        """
        angle = 2 * pi / self.fold
        quaternion = create_quaternion(self.vector, angle)
        return quaternion_rotation(quaternion, coordinate)

    def int_operate(self, coordinate):
        """Returns the rotation symmetry operation on a coordinate
        but returns a tuple of ints.

        This method is useful for seeing how the symmetry operation
        affects the sign of a gaussian functions integral exponents.

        Parameters
        ----------
        coordinate : Tuple[float, float, float]

        Returns
        -------
        : Tuple[int, int, int]

        """
        angle = 2 * pi / self.fold
        quaternion = create_quaternion(self.vector, angle)
        x, y, z = quaternion_rotation(quaternion, coordinate)
        return int(round(x, 1)), int(round(y, 1)), int(round(z, 1))

    @property
    def symmetry_operation(self):
        """Returns the rotation symmetry label.

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
        return 'C_{' + str(round(self.fold, 3)) + '}' + str((i, j, k))
