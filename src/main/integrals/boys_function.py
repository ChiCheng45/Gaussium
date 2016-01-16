from math import exp, gamma
from numba import jit

"""
BoysFunction used in the nuclear-attraction and electron repulsion integrals.

DESCRIPTION
    This class contains the static method, boys_function that approximates the boys function using a summation. The
    methods follows the equation from the Handbook of Computational Chemistry pg 280. The equation is split into two
    parts for either a small or large x. I have included the @jit decoration for this function as its easy to get a
    small speed up here without doing much. Simply install and use anaconda as your python interpreter.

ARGUMENTS

    class BoysFunction:

        Methods
        -------
        @staticmethod
        @jit
        def function(v, x):

            Parameters
            ----------
            v : int
                has different definition depending on which integral is used in
            x : float
                has different definition depending on which integral is used in

            Returns
            -------
            ans : float
                the answer

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
    def calculate(v, x):

        # Approximation of the boys function for small x
        if x <= 25:
            i = 0
            ans = 0
            while 1 > 0:
                seq = (gamma(v + 0.5) / gamma(v + i + 1.5)) * x**i
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
                seq = (gamma(v + 0.5) / gamma(v - i + 1.5)) * x**(-i)
                if seq < 1e-10:
                    break
                ans += seq
                i += 1
            ans *= (1/2) * exp(-x)
            ans = (gamma(v + 0.5) / (2*x**(v + 0.5))) - ans
            return ans

    @staticmethod
    def recursion(v, x, f_v):
        return (exp(-x) + 2 * x * f_v) / (2 * v - 1)
