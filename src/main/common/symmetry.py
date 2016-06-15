import numpy


class Symmetry:

    @classmethod
    def check_sym(cls, basis1, basis2, basis3, basis4):
        r_1 = basis1.coordinates
        r_2 = basis2.coordinates
        r_3 = basis3.coordinates
        r_4 = basis4.coordinates
        l_1 = basis1.integral_exponents
        l_2 = basis2.integral_exponents
        l_3 = basis3.integral_exponents
        l_4 = basis4.integral_exponents

        if cls.even_odd(l_1[0]) * cls.even_odd(l_2[0]) * cls.even_odd(l_3[0]) * cls.even_odd(l_4[0]) == -1:
            if r_1[0] == r_2[0] == r_3[0] == r_4[0]:
                return False
            else:
                return True
        elif cls.even_odd(l_1[1]) * cls.even_odd(l_2[1]) * cls.even_odd(l_3[1]) * cls.even_odd(l_4[1]) == -1:
            if r_1[1] == r_2[1] == r_3[1] == r_4[1]:
                return False
            else:
                return True
        elif cls.even_odd(l_1[2]) * cls.even_odd(l_2[2]) * cls.even_odd(l_3[2]) * cls.even_odd(l_4[2]) == -1:
            if r_1[2] == r_2[2] == r_3[2] == r_4[2]:
                return False
            else:
                return True
        return True

    @staticmethod
    def even_odd(num):
        if num % 2 == 0:
            return 1
        else:
            return -1

    @staticmethod
    def sort(i, j, k, l):
        a = i
        b = j
        c = k
        d = l
        if a > b:
            a, b = b, a
        if c > d:
            c, d = d, c
        if a > c:
            a, c = c, a
            b, d = d, b
        if a == c:
            if b > d:
                a, c = c, a
                b, d = d, b
        return a, b, c, d
