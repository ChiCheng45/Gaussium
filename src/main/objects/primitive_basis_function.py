from math import pi, sqrt
from math import factorial as fac
import numpy as np


class PrimitiveBasis:
    """Primitive Gaussian basis function.

    Attributes
    ----------
    contraction : float
    exponent : float
    coordinates : Tuple[float, float, float]
    integral_exponents : Tuple[int, int, int]
    normalisation_memo : {None, float}
        Stores the normalisation constant once calculated.

    """
    def __init__(self, contraction, exponent, coordinates, integral_exponents):
        self.contraction = contraction
        self.exponent = exponent
        self.coordinates = coordinates
        self.integral_exponents = integral_exponents
        self.normalisation_memo = None

    @property
    def normalisation(self):
        """Calculates normalisation constant and stores in self.normalisation once called.

        Returns
        -------
        self.normalisation_memo : float

        """
        if self.normalisation_memo is None:
            x, y, z = self.integral_exponents
            out1 = ((2 * self.exponent) / pi) ** (3 / 4)
            out2 = (8 * self.exponent) ** (x + y + z) * fac(x) * fac(y) * fac(z)
            out3 = fac(2 * x) * fac(2 * y) * fac(2 * z)
            self.normalisation_memo = out1 * sqrt(out2 / out3)
        return self.normalisation_memo

    def value(self, x, y, z):
        """Returns the value of the function at point x, y, z.

        Parameters
        ----------
        x : float
        y : float
        z : float

        Returns
        -------
        : float

        """
        r_x, r_y, r_z = self.coordinates
        l_x, l_y, l_z = self.integral_exponents
        return self.normalisation * self.contraction * (x - r_x)**l_x * (y - r_y)**l_y * (z - r_z)**l_z \
        * np.exp(-self.exponent * ((x - r_x)**2 + (y - r_y)**2 + (z - r_z)**2))
