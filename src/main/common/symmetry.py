class Symmetry:

    @staticmethod
    def check_sym(gaussian1, gaussian2):
        coordinates1 = gaussian1.coordinates
        coordinates2 = gaussian2.coordinates
        l_1 = gaussian1.integral_exponents
        l_2 = gaussian2.integral_exponents

        x = coordinates1[0] - coordinates2[0]
        y = coordinates1[1] - coordinates2[1]
        z = coordinates1[2] - coordinates2[2]

        if x == 0:
            if (l_1[1] % 2 == 0 and l_2[1] % 2 != 0) or (l_1[1] % 2 != 0 and l_2[1] % 2 == 0):
                return False
            elif (l_1[2] % 2 == 0 and l_2[2] % 2 != 0) or (l_1[2] % 2 != 0 and l_2[2] % 2 == 0):
                return False
        elif y == 0:
            if (l_1[0] % 2 == 0 and l_2[0] % 2 != 0) or (l_1[0] % 2 != 0 and l_2[0] % 2 == 0):
                return False
            elif (l_1[2] % 2 == 0 and l_2[2] % 2 != 0) or (l_1[2] % 2 != 0 and l_2[2] % 2 == 0):
                return False
        elif z == 0:
            if (l_1[0] % 2 == 0 and l_2[0] % 2 != 0) or (l_1[0] % 2 != 0 and l_2[0] % 2 == 0):
                return False
            elif (l_1[1] % 2 == 0 and l_2[1] % 2 != 0) or (l_1[1] % 2 != 0 and l_2[1] % 2 == 0):
                return False
        else:
            return True
