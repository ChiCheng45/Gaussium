from math import sqrt
from src.main.objects import Coordinates

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
    tuple1  input: a tuple of length three
    tuple2  input: a tuple of length three

    dot_product(tuple1, tuple2)
    tuple1  input: a tuple of length three
    tuple2  input: a tuple of length three
    ans     output: the tuple of the dot product of tuple1 and tuple2

    distance(tuple1, tuple2)
    tuple1  input: a tuple of length three
    tuple2  input: a tuple of length three
    r_ab    output: the euclidean distance between the vectors of tuple1 and tuple2

    multiply(x, tuple1)
    tuple1  input: a tuple of length three
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
    def add(r_1, r_2):
        return Coordinates(r_1.x + r_2.x, r_1.y + r_2.y, r_1.z + r_2.z)

    @staticmethod
    def minus(r_1, r_2):
        return Coordinates(r_1.x - r_2.x, r_1.y - r_2.y, r_1.z - r_2.z)

    @staticmethod
    def multiply(a, r_1):
        return Coordinates(a * r_1.x, a * r_1.y, a * r_1.z)

    @staticmethod
    def dot_product(r_1, r_2):
        return r_1.x * r_2.x + r_1.y * r_2.y + r_1.z + r_2.z

    @staticmethod
    def distance(r_1, r_2):
        r_12 = sqrt((r_1.x - r_2.x)**2 + (r_1.y - r_2.y)**2 + (r_1.z - r_2.z)**2)
        return r_12

    @classmethod
    def gaussian(cls, a, r_1, b, r_2):
        i = (a * r_1.x + b * r_2.x) / (a + b)
        j = (a * r_1.y + b * r_2.y) / (a + b)
        k = (a * r_1.z + b * r_2.z) / (a + b)
        return Coordinates(i, j, k)
