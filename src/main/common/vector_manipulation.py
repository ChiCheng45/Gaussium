import math


class VectorManipulation:

    @classmethod
    def squared_distance(cls, tuple1, tuple2):
        tuple3 = tuple([a - b for a, b in zip(tuple1, tuple2)])
        r_ab = math.sqrt(tuple3[0]**2 + tuple3[1]**2 + tuple3[2]**2)
        return r_ab

    @classmethod
    def vector_add(cls, tuple1, tuple2):
        return tuple([a + b for a, b in zip(tuple1, tuple2)])

    @classmethod
    def vector_minus(cls, tuple1, tuple2):
        return tuple([a - b for a, b in zip(tuple1, tuple2)])

    @classmethod
    def vector_multiply(cls, x, tuple1):
        return tuple([a * x for a in tuple1])

    @classmethod
    def vector_gaussian(cls, a, tuple1, b, tuple2):
        tuple3 = cls.vector_multiply(a, tuple1)
        tuple4 = cls.vector_multiply(b, tuple2)
        tuple_final = cls.vector_add(tuple3, tuple4)
        tuple_final = tuple([x / (a + b) for x in tuple_final])
        return tuple_final
