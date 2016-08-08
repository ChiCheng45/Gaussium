from scipy import integrate
import numpy as np


class ExchangeCorrelation:

    def __init__(self, basis_set, exchange_potential, correlation_potential, int_space=10):
        self.basis_set = basis_set
        self.exchange_potential = exchange_potential
        self.correlation_potential = correlation_potential
        self.int_space = int_space
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

        def integrand(x, y, z):
            x_0 = 2 * self.int_space * x - self.int_space
            y_0 = 2 * self.int_space * y - self.int_space
            z_0 = 2 * self.int_space * z - self.int_space
            if i == j:
                return g_i.value(x_0, y_0, z_0)**2 * (self.exchange_potential.calculate(electron_density(x_0, y_0, z_0))
                + self.correlation_potential.calculate(electron_density(x_0, y_0, z_0)))
            else:
                return g_i.value(x_0, y_0, z_0) * (self.exchange_potential.calculate(electron_density(x_0, y_0, z_0))
                + self.correlation_potential.calculate(electron_density(x_0, y_0, z_0))) * g_j.value(x_0, y_0, z_0)

        integral, error = integrate.nquad(integrand, [[0, 1], [0, 1], [0, 1]],
        opts=[{'epsabs': 1e-3, 'epsrel': 0},
              {'epsabs': 1e-3, 'epsrel': 0},
              {'epsabs': 1e-3, 'epsrel': 0}])

        return (2 * self.int_space)**3 * integral
