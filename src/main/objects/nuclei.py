class Nuclei:
    """Nuclei object, containing the nuclei attributes.

    Attributes
    ----------
    element : str
    charge : int
    mass : float
    coordinates : Tuple[float, float, float]

    """
    def __init__(self, element, charge, mass, coordinates):
        self.element = element
        self.charge = charge
        self.mass = mass
        self.coordinates = coordinates
