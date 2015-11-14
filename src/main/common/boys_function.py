from scipy.special import gamma
import numpy as np


class BoysFunction:

    @staticmethod
    def calculate(v, x):
        if x <= 20:
            i = 0
            out2 = 0
            while i > -1:
                out1 = (gamma(v + (1/2)) / gamma(v + i + (3/2))) * x**i
                if out1 < 1e-10:
                    break
                out2 += out1
                i += 1
            out2 *= (1/2) * np.exp(-x)
            return out2
        elif x > 20:
            i = 0
            out2 = 0
            while i > -1:
                out1 = (gamma(v + (1/2)) / gamma(v - i + (3/2))) * x**(-i)
                if out1 < 1e-10:
                    break
                out2 += out1
                i += 1
            out2 *= (1/2) * np.exp(-x)
            out2 = (gamma(v + (1/2)) / (2*x**(v + (1/2)))) - out2
            return out2
