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

        if cls.even_odd(l_1.l) * cls.even_odd(l_2.l) * cls.even_odd(l_3.l) * cls.even_odd(l_4.l) == -1:
            if r_1.x == r_2.x == r_3.x == r_4.x:
                return False
            else:
                return True
        elif cls.even_odd(l_1.m) * cls.even_odd(l_2.m) * cls.even_odd(l_3.m) * cls.even_odd(l_4.m) == -1:
            if r_1.y == r_2.y == r_3.y == r_4.y:
                return False
            else:
                return True
        elif cls.even_odd(l_1.n) * cls.even_odd(l_2.n) * cls.even_odd(l_3.n) * cls.even_odd(l_4.n) == -1:
            if r_1.z == r_2.z == r_3.z == r_4.z:
                return False
            else:
                return True
        else:
            return True

    @staticmethod
    def even_odd(num):
        if num % 2 == 0:
            return 1
        else:
            return -1
