from src.main.kohnsham.exchange import ExchangePotential
import numpy as np


class SlaterExchange(ExchangePotential):

    def __init__(self, alpha):
        super().__init__()
        self.alpha = alpha

    def calculate(self, density):
        """Returns the value of the slater exchange potential for a given density.

        This exchange potential can be used by calling the SLATER keyword.

        Parameters
        ----------
        density : float

        Returns
        -------
        : float

        """
        if density not in self.potential_memo:
            self.potential_memo[density] = -(3/2) * self.alpha * (3/np.pi)**(1/3) * density**(1/3)
        return self.potential_memo[density]
