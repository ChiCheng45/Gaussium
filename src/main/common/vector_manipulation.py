import math


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
        for x in range(0, len(tuple1)):
            ans += tuple1[x] * tuple2[x]
        return ans

    @staticmethod
    def distance(tuple1, tuple2):
        tuple3 = tuple([a - b for a, b in zip(tuple1, tuple2)])
        r_ab = math.sqrt(tuple3[0]**2 + tuple3[1]**2 + tuple3[2]**2)
        return r_ab

    @classmethod
    def gaussian(cls, a, tuple1, b, tuple2):
        tuple3 = cls.multiply(a, tuple1)
        tuple4 = cls.multiply(b, tuple2)
        tuple_final = cls.add(tuple3, tuple4)
        tuple_final = tuple([x / (a + b) for x in tuple_final])
        return tuple_final
