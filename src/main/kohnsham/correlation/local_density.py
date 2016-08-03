from math import log
import numpy as np


def chachiyo_correlation(density):
    wigner_seitz_radius = (4 * np.pi * density / 3)**(-1/3)
    return -0.015545 * log(1 + (20.456255 / wigner_seitz_radius) + (20.456255 / wigner_seitz_radius**2))
