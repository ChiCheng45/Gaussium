from math import exp, gamma


def boys_function(v, x):
    """Computes the boys function used for calculating the two electron and nuclear attraction integrals.

    Parameters
    ----------
    v : {int, float}
    x : {int, float}

    Returns
    -------
    ans : float

    Notes
    -----
    There are also no checks for if the series diverges and causes an infinite loop.

    References
    ----------
    [1] Handbook of Computational Chemistry pg. 280

    """
    # Approximation of the boys function for small x
    if x <= 25:
        i = 0
        ans = 0
        g_v = gamma(v + 0.5)
        while True:
            seq = (g_v / gamma(v + i + 1.5)) * x**i
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
        g_v = gamma(v + 0.5)
        while True:
            seq = (g_v / gamma(v - i + 1.5)) * x**(-i)
            if seq < 1e-10:
                break
            ans += seq
            i += 1
        ans *= (1/2) * exp(-x)
        ans = (g_v / (2*x**(v + 0.5))) - ans
        return ans


def boys_function_recursion(v, x, f_v):
    """Returns the answer to the boys function f_{v - 1}(x) using the answer for the boys function f_{v}(x).

    Parameters
    ----------
    v : {int, float}
    x : {int, float}
    f_v : float

    Returns
    -------
    : float

    References
    ----------
    [1] Handbook of Computational Chemistry pg. 280

    """
    return (exp(-x) + 2 * x * f_v) / (2 * v - 1)
