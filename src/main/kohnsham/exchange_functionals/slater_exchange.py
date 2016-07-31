from scipy import integrate
import numpy as np


class SlaterExchange:

    def __init__(self, basis_set, alpha=1.0):
        self.basis_set = basis_set
        self.alpha = alpha
        self.electron_density_memo = {}
        self.density_matrix_memo = np.matrix([])

    def integrate(self, density_matrix, i, j):

        if density_matrix is not self.density_matrix_memo:
            self.density_matrix_memo = density_matrix
            self.electron_density_memo = {}

        g_i = self.basis_set[i]
        g_j = self.basis_set[j]

        def electron_density(x, y, z):
            if (x, y, z) not in self.electron_density_memo:
                density = 0
                for a, basis_a in enumerate(self.basis_set):
                    for b, basis_b in enumerate(self.basis_set):
                        if a == b:
                            density += density_matrix.item(a, b) * basis_a.value(x, y, z)**2
                        elif a <= b:
                            density += 2 * basis_a.value(x, y, z) * density_matrix.item(a, b) * basis_b.value(x, y, z)
                self.electron_density_memo[(x, y, z)] = density
            return self.electron_density_memo[(x, y, z)]

        def integrand_ii(x, y, z):
            return (electron_density(x, y, z))**(1/3) * g_i.value(x, y, z)**2

        def integrand_ij(x, y, z):
            return g_i.value(x, y, z) * (electron_density(x, y, z))**(1/3) * g_j.value(x, y, z)

        if i == j:
            integral, error = integrate.nquad(integrand_ii, [[-6, 6], [-6, 6], [-6, 6]],
            opts=[{'epsabs': 0, 'epsrel': 1e-3}, {'epsabs': 0, 'epsrel': 1e-3}, {'epsabs': 0, 'epsrel': 1e-3}])
        else:
            integral, error = integrate.nquad(integrand_ij, [[-6, 6], [-6, 6], [-6, 6]],
            opts=[{'epsabs': 0, 'epsrel': 1e-3}, {'epsabs': 0, 'epsrel': 1e-3}, {'epsabs': 0, 'epsrel': 1e-3}])

        return - self.alpha * (3/2) * (3/np.pi)**(1/3) * integral
