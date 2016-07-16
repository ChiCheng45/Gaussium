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
