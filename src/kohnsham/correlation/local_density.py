from src.kohnsham.correlation import CorrelationPotential
import numpy as np
from numpy import vectorize


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
        self.a = a
        self.x_0 = x_0
        self.b = b
        self.c = c

    def potential(self, density):
        vfunc = vectorize(self._potential)
        return vfunc(density)

    def _potential(self, density):
        """Returns the value of the Vosko-Wilk-Nusair correlation potential for a given density and fit parameters.

        Parameters
        ----------
        density : float

        Returns
        -------
        : float

        """
        if density == 0:
            return 0
        else:
            a = self.a
            x_0 = self.x_0
            b = self.b
            c = self.c
            x = self.wigner_seitz_radius(density)**(1/2)
            x_x = x**2 + b * x + c
            x_x_0 = x_0**2 + b * x_0 + c
            q = (4 * c - b**2)**(1/2)

            return a * (np.log(x**2 / x_x) + (2 * b / q) * np.arctan(q / (2 * x + b))
            - (b * x_0 / x_x_0) * (np.log((x - x_0)**2 / x_x) + (2 * (b + 2 * x_0) / q) * np.arctan(q / (2 * x + b))))
