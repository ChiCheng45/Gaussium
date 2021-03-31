import numpy as np
from src.kohnsham.exchange import ExchangePotential


class SlaterExchange(ExchangePotential):
    """Slater exchange potential for a given value alpha.

    Attributes
    ----------
    alpha : float

    """
    def __init__(self, alpha):
        self.alpha = alpha

    def energy(self, density):
        """Returns the value of the slater exchange energy per electron
        for a given density.

        This exchange potential can be used by calling the 'S' or
        'XA' keyword.

        Parameters
        ----------
        density : float

        Returns
        -------
        : float

        """
        return -(9/8) * self.alpha * (3/np.pi)**(1/3) * density**(1/3)

    def potential(self, density):
        """Returns the value of the slater exchange potential
        for a given density.

        This exchange potential can be used by calling the 'S' or
        'XA' keyword.

        Parameters
        ----------
        density : float

        Returns
        -------
        : float

        """
        return -(3/2) * self.alpha * (3/np.pi)**(1/3) * density**(1/3)
