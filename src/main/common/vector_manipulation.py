from math import sqrt


"""
NAME
    Vector

SYNOPSIS
    def add(tuple1, tuple2)
    def minus(tuple1, tuple2)
    def multiply(x, tuple1)
    def dot_product(tuple1, tuple2)
    def distance(tuple1, tuple2)
    def gaussian(x, tuple1, y, tuple2)
    def tuple tuple1, tuple2
    double x, y

DESCRIPTION
    The Vector class acts as a container for all the basic tuple manipulation. As the class stores the coordinates as
    tuples we use these methods to carry out all the important mathematics. For some reason, python is far faster having
    these methods then to use the scipy or numpy equivalents.
    All the below methods are obvious except gaussian method which is a method which returns the coordinates of a
    gaussian function from the product of two other gaussian function with coordinates tuple1, tuple2 and with exponents
    x and y.

ARGUMENTS
    def add(tuple1, tuple2), minus(tuple1, tuple2)
    tuple1  input: a tuple of length three
    tuple2  input: a tuple of length three

    def dot_product(tuple1, tuple2)
    tuple1  input: a tuple of length three
    tuple2  input: a tuple of length three
    ans     output: the tuple of the dot product of tuple1 and tuple2

    def distance(tuple1, tuple2)
    tuple1  input: a tuple of length three
    tuple2  input: a tuple of length three
    r_ab    output: the euclidean distance between the vectors of tuple1 and tuple2

    def multiply(x, tuple1)
    tuple1  input: a tuple of length three
    x       input: x any number

    def gaussian(x, tuple1, y, tuple2)
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
        return r_1[0] + r_2[0], r_1[1] + r_2[1], r_1[2] + r_2[2]

    @staticmethod
    def minus(r_1, r_2):
        return r_1[0] - r_2[0], r_1[1] - r_2[1], r_1[2] - r_2[2]

    @staticmethod
    def multiply(a, r_1):
        return a * r_1[0], a * r_1[1], a * r_1[2]

    @staticmethod
    def dot_product(r_1, r_2):
        return r_1[0] * r_2[0] + r_1[1] * r_2[1] + r_1[2] + r_2[2]

    @staticmethod
    def distance(r_1, r_2):
        r_12 = sqrt((r_1[0] - r_2[0])**2 + (r_1[1] - r_2[1])**2 + (r_1[2] - r_2[2])**2)
        return r_12

    @staticmethod
    def gaussian_product(a, r_1, b, r_2):
        i = (a * r_1[0] + b * r_2[0]) / (a + b)
        j = (a * r_1[1] + b * r_2[1]) / (a + b)
        k = (a * r_1[2] + b * r_2[2]) / (a + b)
        return i, j, k
