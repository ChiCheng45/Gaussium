from math import pi, sqrt
from scipy.special import factorial2
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
        """Calculates normalisation constant for the primitive gaussian
        and stores in self.normalisation once called.

        Returns
        -------
        self.normalisation_memo : float

        """
        if self.normalisation_memo is None:
            l, m, n = self.integral_exponents
            out1 = factorial2(2 * l - 1) * factorial2(2 * m - 1) * factorial2(2 * n - 1)
            out2 = (pi / (2 * self.exponent))**(3/2)
            out3 = (4 * self.exponent)**(l + m + n)
            self.normalisation_memo = 1 / sqrt((out1 * out2) / out3)
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
        * np.exp(- self.exponent * ((x - r_x)**2 + (y - r_y)**2 + (z - r_z)**2))
