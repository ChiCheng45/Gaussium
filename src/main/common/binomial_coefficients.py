import numpy as np


class BinomialCoefficientsFunction:

    @classmethod
    def calculate_coefficient(cls, j, l, m, a, b):
        coefficient = 0
        for k in range(max(0, j - m), min(j, l) + 1):
            coefficient += cls.combination(l, k) * cls.combination(m, j - k) * a**(l - k) * b**(m + k - j)
        return coefficient

    @staticmethod
    def combination(n, k):
        if k <= n:
            combination = np.math.factorial(n) / (np.math.factorial(k) * np.math.factorial(n - k))
            return combination
        else:
            return 0
