from scipy.special import factorial2
from math import sqrt, pi
import itertools


class Basis:
    """Contracted Gaussian basis function.

    Attributes
    ----------
    primitive_gaussian_array : List[PrimitiveBasis]
    coordinates : Tuple[float, float, float]
    integral_exponents : Tuple[int, int, int]

    """
    def __init__(self, primitive_gaussian_array, coordinates, integral_exponents):
        self.primitive_gaussian_array = primitive_gaussian_array
        self.coordinates = coordinates
        self.integral_exponents = integral_exponents
        self.normalisation_memo = None

    @property
    def normalisation(self):
        """Calculates normalisation constant for the basis set and stores in self.normalisation once called.

        Returns
        -------
        self.normalisation_memo : float

        """
        if self.normalisation_memo is None:
            l, m, n = self.integral_exponents

            if len(self.primitive_gaussian_array) == 1:
                self.normalisation_memo = 1
            else:
                ans = 0.0
                for primitive_a, primitive_b in itertools.product(self.primitive_gaussian_array, repeat=2):
                    a_1 = primitive_a.exponent
                    a_2 = primitive_b.exponent
                    c_1 = primitive_a.contraction
                    c_2 = primitive_b.contraction
                    n_1 = primitive_a.normalisation
                    n_2 = primitive_b.normalisation

                    out1 = factorial2(2 * l - 1) * factorial2(2 * m - 1) * factorial2(2 * n - 1)
                    out2 = (pi / (a_1 + a_2)) ** (3 / 2)
                    out3 = (2 * (a_1 + a_2)) ** (l + m + n)
                    ans += (c_1 * c_2 * n_1 * n_2 * out1 * out2) / out3
                self.normalisation_memo = 1 / sqrt(ans)

        return self.normalisation_memo

    def value(self, x, y, z):
        """Returns the value at point x, y, z.

        Parameters
        ----------
        x : float
        y : float
        z : float

        Returns
        -------
        ans : float
        """
        ans = 0
        for primitive in self.primitive_gaussian_array:
            ans += primitive.value(x, y, z)
        return self.normalisation * ans