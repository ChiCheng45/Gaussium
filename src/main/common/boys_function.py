from scipy.special import gamma
import numpy as np

"""
NAME
    BoysFunction

SYNOPSIS
    boys_function(v, x)
    int v
    double x

DESCRIPTION
    This class contains the static method, boys_function that approximates the boys function using a summation. The
    methods follows the equation from the Handbook of Computational Chemistry pg 280. The equation is split into two
    parts for either a small or large x.

ARGUMENTS
    boys_function(v, x)
    v       input: has different definition depending on which integral is used in.
    x       input: has different definition depending on which integral is used in.
    ans     output: the answer.

SEE ALSO
    nuclear_attraction_integral.py
    two_electron_repulsion_integral.py

DIAGNOSTICS
    Potential for the while loops to go on till infinity if the series diverges.
"""


class BoysFunction:

    @staticmethod
    def function(v, x):
        # Approximation of the boys function for small x
        if x <= 20:
            i = 0
            ans = 0
            while i > -1:
                seq = (gamma(v + (1/2)) / gamma(v + i + (3/2))) * x**i
                if seq < 1e-10:
                    break
                ans += seq
                i += 1
            ans *= (1/2) * np.exp(-x)
            return ans

        # Approximation of the boys function for large x
        elif x > 20:
            i = 0
            ans = 0
            while i > -1:
                seq = (gamma(v + (1/2)) / gamma(v - i + (3/2))) * x**(-i)
                if seq < 1e-10:
                    break
                ans += seq
                i += 1
            ans *= (1/2) * np.exp(-x)
            ans = (gamma(v + (1/2)) / (2*x**(v + (1/2)))) - ans
            return ans
