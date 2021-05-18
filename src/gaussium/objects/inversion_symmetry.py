class InversionSymmetry:
    """Inversion symmetry operator.

    Attributes
    ----------
    symmetry_operation : str
        The symmetry operation label.

    """
    def __init__(self):
        self.symmetry_operation = 'i'

    def operate(self, coordinate):
        """Returns the inversion symmetry operation on a coordinate
        through the origin.

        Parameters
        ----------
        coordinate : Tuple[float, float, float]

        Returns
        -------
        : Tuple[float, float, float]

        """
        x, y, z = coordinate
        return -x, -y, -z

    def int_operate(self, coordinate):
        """Returns the inversion symmetry operation on a coordinate
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
        x, y, z = coordinate
        return -int(round(x, 1)), -int(round(y, 1)), -int(round(z, 1))
