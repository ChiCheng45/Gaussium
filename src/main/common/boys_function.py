from math import exp, gamma
from numba import jit

"""
BoysFunction

SYNOPSIS
    def boys_function(v, x)
    int v
    double x

DESCRIPTION
    This class contains the static method, boys_function that approximates the boys function using a summation. The
    methods follows the equation from the Handbook of Computational Chemistry pg 280. The equation is split into two
    parts for either a small or large x. I have included the @jit decoration for this function as its easy to get a
    small speed up here without doing much. Simply install and use anaconda as your python interpreter.

ARGUMENTS
    def boys_function(v, x)
    v       input: has different definition depending on which integral is used in.
    x       input: has different definition depending on which integral is used in.
    ans     output: the answer.

SEE ALSO
    nuclear_attraction_integral.py
    cook_integral.py
    http://numba.pydata.org/

DIAGNOSTICS
    Potential for the while loops to go on to infinity if the series diverges.
"""


class BoysFunction:

    @staticmethod
    @jit
    def function(v, x):
        # Approximation of the boys function for small x
        if x <= 25:
            i = 0
            ans = 0
            while 1 > 0:
                seq = (gamma(v + (1/2)) / gamma(v + i + (3/2))) * x**i
                if seq < 1e-10:
                    break
                ans += seq
                i += 1
            ans *= (1/2) * exp(-x)
            return ans

        # Approximation of the boys function for large x
        elif x > 25:
            i = 0
            ans = 0
            while 1 > 0:
                seq = (gamma(v + (1/2)) / gamma(v - i + (3/2))) * x**(-i)
                if seq < 1e-10:
                    break
                ans += seq
                i += 1
            ans *= (1/2) * exp(-x)
            ans = (gamma(v + (1/2)) / (2*x**(v + (1/2)))) - ans
            return ans
