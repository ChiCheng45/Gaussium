from math import factorial as fac

"""
Binomial class containing methods for calculating the coefficients of (x+a)**l * (x+b)**m type equations.

DESCRIPTION
    This class is a container for the methods used for calculating the binomial coefficients. The combination static
    method simply returns the combination from mathematics. While the class method calculate_coefficients calculates the
    coefficient for x**j from the expansion of

        (x+a)**l * (x+b)**m  -- (1)

    The method coefficients uses the method combination to shorten the code a little. The formula of coefficients can be
    found in the Handbook of Computational Chemistry pg 237. These functions are important for the molecular integrals
    of p-type functions and higher as we need to multiply eqn.1 type functions.

ARGUMENTS

    class Binomial:

        Methods
        -------
        @staticmethod
        def combination(n, k):

            Parameters
            ----------
            n, k : int
                simply the definitions from mathematics

            Returns
            -------
            combination : int
                the output of the combination from mathematics

        @staticmethod
        def coefficient(cls, j, l_1, l_2, a, b):

            Parameters
            ----------
            j : int
                determines the power of x you want to calculate the coefficient of, e.g.
                coefficient(cls, 2, l_1, l_2, a, b) will calculate the coefficient for x**2
            l_1, l_2 : int
                these are the powers from eqn.1. When used in the molecular integral these will end up being the l or m
                or n exponents from two different gaussian functions
            a, b : float
                a and b values of eqn.1. When used in the molecular integrals this will become the PA_x, PB_x or PA_y,
                PB_y or PA_z, PB_z values where the vector P is formed from multiplying two gaussian functions centered
                on coordinates A and B

            Returns
            -------
            coefficient : float
                the coefficient of x**j from multiplying out eqn.1

SEE ALSO
    nuclear_attraction_integral.py
    orbital_overlap_integral.py
    cook_integral.py
    https://en.wikipedia.org/wiki/Combination

DIAGNOSTICS
    None
"""


class Binomial:

    @staticmethod
    def combination(n, k):
        if k <= n:
            ans = fac(n) / (fac(k) * fac(n - k))
            return ans
        else:
            return 0

    @classmethod
    def coefficient(cls, j, l_1, l_2, a, b):
        ans = 0
        for k in range(max(0, j - l_2), min(j, l_1) + 1):
            ans += cls.combination(l_1, k) * cls.combination(l_2, j - k) * a**(l_1 - k) * b**(l_2 + k - j)
        return ans
