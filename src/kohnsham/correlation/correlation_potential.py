import numpy as np


class CorrelationPotential:

    def wigner_seitz_radius(self, density):
        """Returns the wigner sietz radius for a given density.

        Parameters
        ----------
        density : float

        Returns
        -------
        : ans

        """
        return (3 / (4 * np.pi * density))**(1/3)

    def calculate(self, density):
        """Returns 0.0 for cases where no correlation potential was selected.

        Parameters
        ----------
        density : float

        Returns
        -------
        : float

        """
        return 0.0
