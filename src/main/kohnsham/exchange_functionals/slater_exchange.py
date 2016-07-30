from scipy import integrate
import numpy as np


class SlaterExchange:

    def __init__(self, basis_set):
        self.basis_set = basis_set

    def integrate(self, density_matrix, i, j):

        g_i = self.basis_set[i]
        g_j = self.basis_set[j]

        def electron_density(x, y, z):
            density = 0
            for a, basis_a in enumerate(self.basis_set):
                for b, basis_b in enumerate(self.basis_set):
                    density += basis_a.value(x, y, z) * density_matrix.item(a, b) * basis_b.value(x, y, z)
            return density

        def integrand(x, y, z):
            return g_i.value(x, y, z) * (electron_density(x, y, z))**(1/3) * g_j.value(x, y, z)

        integral, error = integrate.nquad(integrand, [[-3, 3], [-3, 3], [-3, 3]], opts=[{'epsabs': 0, 'epsrel': 1e-3},
        {'epsabs': 0, 'epsrel': 1e-3}, {'epsabs': 0, 'epsrel': 1e-3}])
        return - (3/2) * (3/np.pi)**(1/3) * integral
