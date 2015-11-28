from math import factorial as fac

"""
NAME
    BinomialCoefficientsFunction

SYNOPSIS
    combination(n, k)
    int n, k, combination

    calculate_coefficient(cls, j, l_1, l_2, a, b)
    int j, l_1, l_2
    float a, b, coefficient

DESCRIPTION
    This class is a container for the functions used for calculating the binomial coefficients. The combination static
    method simply returns the combination from mathematics. While the class method calculate_coefficients calculates the
    coefficient for x**j from the expansion of

        (x+a)**l * (x+b)**m  -- (1)

    The method calculate_coefficients uses the method combination to shorten the code a little. The formula can be found
    of calculate_coefficients in the Handbook of Computational Chemistry pg 237. These functions are important for the
    molecular integrals for p-type and higher as we need to multiply the polynomials eqn.1.

ARGUMENTS
    combination(n, k)
    n, k        Input:  simply the definitions from mathematics
    combination Output: the output of the combination from mathematics

    calculate_coefficient(cls, j, l, m, a, b)
    j           Input:  determine the power of x you want to calculate the coefficient for, e.g.
                        calculate_coefficient(cls, 2, l, m, a, b) will calculated the coefficient for x**2
    l_1, l_2    Input:  these are the powers of the polynomials from eqn.1. When used in the molecular integral
                        these will end up being the l or m or n exponents from two gaussian functions.
    a, b        Input:  a and b values of eqn.1. When used in the molecular integrals this will become the PA_x, PB_x or
                        PA_y, PB_y or PA_z, PB_z values where the vector P is formed from multiplying two gaussian
                        functions centered on coordinates A and B.
    coefficient Output: the coefficient of x**j for eqn.1.

SEE ALSO
    nuclear_attraction_integral.py
    orbital_overlap_integral.py
    two_electron_repulsion_integral.py
    https://en.wikipedia.org/wiki/Combination

DIAGNOSTICS
    None
"""


class Binomial:

    @staticmethod
    def combination(n, k):
        if k <= n:
            combination = fac(n) / (fac(k) * fac(n - k))
            return combination
        else:
            return 0

    @staticmethod
    def calculate_coefficient(j, l_1, l_2, a, b):
        coefficient = 0
        for k in range(max(0, j - l_2), min(j, l_1) + 1):
            coefficient += Binomial.combination(l_1, k) * Binomial.combination(l_2, j - k) * a**(l_1 - k) * b**(l_2 + k - j)
        return coefficient
