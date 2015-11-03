from math import factorial, sqrt
from scipy.special import erf
import numpy as np


class GammaFunction:

    @staticmethod
    def incomplete_gamma_function(v, u):
        if u == 0:
            out = 1
        else:
            div = factorial(2*v) / (2*factorial(v))
            left = ((sqrt(np.pi) * erf(sqrt(u))) / (4**v * u**(v + (1/2))))
            right = 0
            for k in range(0, v):
                right += factorial(v - k) / (4**k * factorial(2*v - 2*k) * u**(k + 1))
            right *= np.exp(-u)
            out = div * (left - right)
        return out
