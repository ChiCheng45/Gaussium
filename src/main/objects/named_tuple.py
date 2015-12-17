from collections import namedtuple


class Coordinates(namedtuple('Coordinates', ('x', 'y', 'z'))):
    pass


class IntegralExponents(namedtuple('IntegralExponents', ('l', 'm', 'n'))):
    pass
