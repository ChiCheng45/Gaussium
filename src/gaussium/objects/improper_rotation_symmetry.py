from math import pi
from gaussium.common import create_householder_matrix
from gaussium.common import create_quaternion
from gaussium.common import householder_matrix_reflection
from gaussium.common import quaternion_rotation


class ImproperRotationSymmetry:
    """Improper rotation symmetry operator.

    This class can only create axis of improper rotation that goes
    through the origin. As the reflection plane of the improper rotation
    will go through the origin, it is therefore equivalent to rotation
    and inversion. This operator will instead rotate and invert as the
    code is simpler and does not need the creation of the householder
    matrix used in reflection operations.

    Attributes
    ----------
    fold : float
        The fraction of 2 * pi angle a coordinate is rotated.
    vector : Tuple[float, float, float]
        The axis of improper rotation.

    """
    def __init__(self, fold, vector):
        self.fold = fold
        self.vector = vector

    def operate(self, coordinate):
        """Rotates a point around the axis of rotation for a angle of
        2 * pi / self.fold and then inverts.

        Parameters
        ----------
        coordinate : Tuple[float, float, float]

        Returns
        -------
        : Tuple[float, float, float]

        """
        angle = 2 * pi / self.fold
        quaternion = create_quaternion(self.vector, angle)
        householder_matrix = create_householder_matrix(self.vector)
        x, y, z = householder_matrix_reflection(quaternion_rotation(quaternion, coordinate), householder_matrix)
        return x, y, z

    def int_operate(self, coordinate):
        """Returns the rotation symmetry operation on a coordinate but
        returns a tuple of ints.

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
        householder_matrix = create_householder_matrix(self.vector)
        x, y, z = householder_matrix_reflection(quaternion_rotation(quaternion, coordinate), householder_matrix)
        return -int(round(x, 1)), -int(round(y, 1)), -int(round(z, 1))

    @property
    def symmetry_operation(self):
        """Returns the improper rotation symmetry label.

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
        return 'S_{' + str(round(self.fold, 3)) + '}' + str((i, j, k))
