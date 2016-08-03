import numpy as np


def slater_exchange(density):
    return -(3/2) * (3/np.pi)**(1/3) * density**(1/3)
