class ExchangePotential:

    def __init__(self):
        self.potential_memo = {}

    def calculate(self, density):
        """Returns 0.0 for cases where no exchange potential was selected.

        Parameters
        ----------
        density : float

        Returns
        -------
        : float

        """
        return 0.0
