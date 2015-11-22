from math import sqrt

"""
NAME
    Vector

SYNOPSIS
    add(tuple1, tuple2)
    minus(tuple1, tuple2)
    multiply(x, tuple1)
    dot_product(tuple1, tuple2)
    distance(tuple1, tuple2)
    gaussian(cls, x, tuple1, y, tuple2)
    tuple tuple1, tuple2
    double x, y

DESCRIPTION
    The Vector class acts as a container for all the basic tuple manipulation. As the class stores the coordinates as
    tuples we use these methods to carry out all the important mathematics. For some reason, python is far faster having
    these methods then to use the scipy or numpy equivalents.
    All the below methods are obvious except gaussian method which is a method which returns the coordinates of a
    gaussian function from the product of two other gaussian function with coordinates tuple1, tuple2 and with exponents
    x and y.

ARGUMENTS
    add(tuple1, tuple2), minus(tuple1, tuple2)
    tuple1  input: a tuple of any length
    tuple2  input: a tuple of any length

    dot_product(tuple1, tuple2)
    tuple1  input: a tuple of any length
    tuple2  input: a tuple of any length
    ans     output: the tuple of the dot product of tuple1 and tuple2

    distance(tuple1, tuple2)
    tuple1  input: a tuple of any length
    tuple2  input: a tuple of any length
    r_ab    output: the euclidean distance between the vectors of tuple1 and tuple2

    multiply(x, tuple1)
    tuple1  input: a tuple of any length
    x       input: x any number

    gaussian(cls, x, tuple1, y, tuple2)
    tuple1  input: the coordinates of gaussian1
    tuple2  input: the coordinates of gaussian2
    x       input: the exponent of gaussian1
    y       input: the exponent of gaussian2
    ans     output: a tuple of the coordinates of a gaussian made from the product of gaussian1 and gaussian2

SEE ALSO
    nuclear_attraction_integral.py
    orbital_overlap_integral.py
    two_electron_repulsion_integral.py

DIAGNOSTICS
    None

"""


class Vector:

    @staticmethod
    def add(tuple1, tuple2):
        return tuple([a + b for a, b in zip(tuple1, tuple2)])

    @staticmethod
    def minus(tuple1, tuple2):
        return tuple([a - b for a, b in zip(tuple1, tuple2)])

    @staticmethod
    def multiply(x, tuple1):
        return tuple([a * x for a in tuple1])

    @staticmethod
    def dot_product(tuple1, tuple2):
        ans = 0
        for x in range(len(tuple1)):
            ans += tuple1[x] * tuple2[x]
        return ans

    @staticmethod
    def distance(tuple1, tuple2):
        tuple3 = tuple([a - b for a, b in zip(tuple1, tuple2)])
        r_ab = sqrt(tuple3[0]**2 + tuple3[1]**2 + tuple3[2]**2)
        return r_ab

    @classmethod
    def gaussian(cls, x, tuple1, y, tuple2):
        tuple3 = cls.multiply(x, tuple1)
        tuple4 = cls.multiply(y, tuple2)
        ans = cls.add(tuple3, tuple4)
        ans = tuple([a / (x + y) for a in ans])
        return ans
