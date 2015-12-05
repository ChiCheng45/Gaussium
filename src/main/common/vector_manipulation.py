from math import sqrt
from numba import jit

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
    cook_integral.py

DIAGNOSTICS
    None

"""


class Vector:

    @staticmethod
    def add(a, b):
        return a[0] + b[0], a[1] + b[1], a[2] + b[2]

    @staticmethod
    def minus(a, b):
        return a[0] - b[0], a[1] - b[1], a[2] - b[2]

    @staticmethod
    def multiply(x, a):
        return a[0] * x, a[1] * x, a[2] * x

    @staticmethod
    def dot_product(a, b):
        return a[0] * b[0] + a[1] * b[1] + a[2] + b[2]

    @staticmethod
    def distance(a, b):
        r_ab = sqrt((a[0] - b[0])**2 + (a[1] - b[1])**2 + (a[2] - b[2])**2)
        return r_ab

    @staticmethod
    def gaussian(x, a, y, b):
        c = Vector.multiply(x, a)
        d = Vector.multiply(y, b)
        ans = Vector.add(c, d)
        return ans[0] / (x + y), ans[1] / (x + y), ans[2] / (x + y)
