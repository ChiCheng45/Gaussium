import itertools


class Indices:

    def __init__(self, occupied, unoccupied):
        self.occupied_indices = range(occupied)
        self.unoccupied_indices = range(occupied, occupied + unoccupied)

    def singles(self):
        for i in self.occupied_indices:
            for a in self.unoccupied_indices:
                yield i, a

    def doubles(self):
        for i, j in itertools.permutations(self.occupied_indices, 2):
            for a, b in itertools.permutations(self.unoccupied_indices, 2):
                yield i, j, a, b

    def triples(self):
        for i, j, k in itertools.permutations(self.occupied_indices, 3):
            for a, b, c in itertools.permutations(self.unoccupied_indices, 3):
                yield i, j, k, a, b, c

    def restricted_doubles(self):
        for i, j in itertools.product(self.occupied_indices, repeat=2):
            for a, b in itertools.product(self.unoccupied_indices, repeat=2):
                yield i, j, a, b
