import math


class VectorManipulation:

    @staticmethod
    def squared_distance(tuple1, tuple2):
        tuple3 = tuple([a - b for a, b in zip(tuple1, tuple2)])
        r_ab = math.sqrt(tuple3[0]**2 + tuple3[1]**2 + tuple3[2]**2)
        return r_ab

    @staticmethod
    def vector_add(tuple1, tuple2):
        return tuple([a + b for a, b in zip(tuple1, tuple2)])

    @staticmethod
    def vector_minus(tuple1, tuple2):
        return tuple([a - b for a, b in zip(tuple1, tuple2)])

    @staticmethod
    def vector_multiply(x, tuple1):
        return tuple([a * x for a in tuple1])

    @staticmethod
    def dot_product(tuple1, tuple2):
        ans = 0
        for x in range(0, len(tuple1)):
            ans += tuple1[x] * tuple2[x]
        return ans

    @classmethod
    def vector_gaussian(cls, a, tuple1, b, tuple2):
        tuple3 = cls.vector_multiply(a, tuple1)
        tuple4 = cls.vector_multiply(b, tuple2)
        tuple_final = cls.vector_add(tuple3, tuple4)
        tuple_final = tuple([x / (a + b) for x in tuple_final])
        return tuple_final
