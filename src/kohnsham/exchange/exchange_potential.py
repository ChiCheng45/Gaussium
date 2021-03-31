class ExchangePotential:

    def energy(self, density):
        """Returns 0.0 for cases where no exchange energy per electron
        was selected.

        Parameters
        ----------
        density : float

        Returns
        -------
        : float

        """
        return 0.0

    def potential(self, density):
        """Returns 0.0 for cases where no exchange potential
        was selected.

        Parameters
        ----------
        density : float

        Returns
        -------
        : float

        """
        return 0.0
