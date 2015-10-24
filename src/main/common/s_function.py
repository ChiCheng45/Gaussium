from scipy.misc import factorial2


class SFunction:

    def __init__(self, binomial_coefficient):
        self.binomial_coefficient = binomial_coefficient

    def calc(self, l_1, l_2, a, b, gamma):
        s = 0
        for j in range(0, int((l_1 + l_2) / 2) + 1):
            s += self.binomial_coefficient().calculate_coefficient(2*j, l_1, l_2, a, b) * (factorial2(2*j - 1) / (2*gamma)**j)
        return s
