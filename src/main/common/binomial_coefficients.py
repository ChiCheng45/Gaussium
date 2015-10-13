import numpy as np


class BinomialCoefficientsFunction:

    def calculate_coefficient(self, j, l, m, a, b):
        coefficient = 0
        for k in range(max(0, j - m), min(j, l) + 1):
            coefficient += self.combination(l, k) * self.combination(m, j - k) * a**(l - k) * b**(m + k - j)
        return coefficient

    def combination(self, n, k):
        if k <= n:
            combination = np.math.factorial(n) / np.math.factorial(k) * np.math.factorial(n - k)
            return combination
        else:
            return 0
