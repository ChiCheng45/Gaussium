from math import exp, gamma
from numba import jit


@jit
def boys_function(v, x):

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


def boys_function_recursion(v, x, f_v):
    return (exp(-x) + 2 * x * f_v) / (2 * v - 1)
