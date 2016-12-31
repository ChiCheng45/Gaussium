import itertools


class Indices:
    """Convenience functions that yield indices for excitation of a electron from an occupied to unoccupied orbital.

    Correlated methods such as MP2 and CCSD can inherit from this class and can use the methods in the following
    manner,

        def mp2_initial_guess(self):
            t = {}
            for i, a in self.singles():
                t[i, a] = 0
            for i, j, a, b in self.doubles():
                t[i, j, a, b] = self.integrals.item(i, j, a, b) / self.denominator[i, j, a, b]
            return t

    this example shows mp2 guess singles and doubles amplitudes calculated for CCSD calculations.

    Attributes
    ----------
    occupied_indices : int
        number of occupied orbitals
    unoccupied_indices : int
        number of unoccupied orbitals

    """
    def __init__(self, occupied, unoccupied):
        self.occupied_indices = range(occupied)
        self.unoccupied_indices = range(occupied, occupied + unoccupied)

    def singles(self):
        """Lazy evaluation of singles indices.

        Returns
        -------
        i, a : Tuple[int]

        """
        for i in self.occupied_indices:
            for a in self.unoccupied_indices:
                yield i, a

    def doubles(self):
        """Lazy evaluation of doubles indices.

        Returns
        -------
        i, j, a, b : Tuple[int]

        """
        for i, j in itertools.permutations(self.occupied_indices, 2):
            for a, b in itertools.permutations(self.unoccupied_indices, 2):
                yield i, j, a, b

    def triples(self):
        """Lazy evaluation of triples indices.

        Returns
        -------
        i, j, k, a, b, c : Tuple[int]

        """
        for i, j, k in itertools.permutations(self.occupied_indices, 3):
            for a, b, c in itertools.permutations(self.unoccupied_indices, 3):
                yield i, j, k, a, b, c

    def restricted_doubles(self):
        """Lazy evaluation of restricted doubles indices.

        Returns
        -------
        i, j, a, b : Tuple[int]

        """
        for i, j in itertools.product(self.occupied_indices, repeat=2):
            for a, b in itertools.product(self.unoccupied_indices, repeat=2):
                yield i, j, a, b
