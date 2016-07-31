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
        self.value_memo = {}

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
        if (x, y, z) not in self.value_memo:
            ans = 0
            for primitive in self.primitive_gaussian_array:
                ans += primitive.value(x, y, z)
            self.value_memo[(x, y, z)] = ans
        return self.value_memo[(x, y, z)]
