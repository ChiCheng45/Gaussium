from src.main.kohnsham.correlation import CorrelationPotential
from math import log, atan


class VoskoWilkNusair(CorrelationPotential):
    """Vosko-Wilk-Nusair correlation potential.

    Correlation potential object can be initialized with any fitting parameters, VWN3 and VWN5 examples below.

    Attributes
    ----------
    a : float
        VWN3 paramagnetic = 0.0621814
        VWN5 paramagnetic = 0.0621814
        VWN3 ferromagnetic = 0.0310907
        VWN5 ferromagnetic = 0.0310907
    x_0 : float
        VWN3 paramagnetic = -0.409286
        VWN5 paramagnetic = -0.10498
        VWN3 ferromagnetic = -0.743294
        VWN5 ferromagnetic = -0.32500
    b : float
        VWN3 paramagnetic = 13.0720
        VWN5 paramagnetic = 3.72744
        VWN3 ferromagnetic = 20.1231
        VWN5 ferromagnetic = 7.06042
    c : float
        VWN3 paramagnetic = 42.7198
        VWN5 paramagnetic = 12.9352
        VWN3 ferromagnetic = 101.578
        VWN5 ferromagnetic = 18.0578

    References
    ----------
    S. H. Vosko, L. Wilk, M. Nusair, Canadian Journal of Physics, 1980, 58(8): 1200-1211, 10.1139/p80-159

    """
    def __init__(self, a, x_0, b, c):
        super().__init__()
        self.a = a
        self.x_0 = x_0
        self.b = b
        self.c = c

    def calculate(self, density):
        """Returns the value of the Vosko-Wilk-Nusair correlation potential for a given density and fit parameters.

        Parameters
        ----------
        density : float

        Returns
        -------
        : float

        """
        if density not in self.potential_memo:
            a = self.a
            x_0 = self.x_0
            b = self.b
            c = self.c
            x = self.wigner_seitz_radius(density)**(1/2)
            x_x = x**2 + b * x + c
            x_x_0 = x_0**2 + b * x_0 + c
            q = (4 * c - b**2)**(1/2)

            self.potential_memo[density] = a * (log(x**2 / x_x) + (2 * b / q) * atan(q / (2 * x + b))
            - (b * x_0 / x_x_0) * (log((x - x_0)**2 / x_x) + (2 * (b + 2 * x_0) / q) * atan(q / (2 * x + b)))) \
            - (a / 3) * ((1 + x_0 * x) / (1 + x_0 * x + b * x**2 + c * x**3))

        return self.potential_memo[density]
