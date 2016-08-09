from scipy import integrate
import numpy as np


class ExchangeCorrelation:

    def __init__(self, basis_set, exchange_potential, correlation_potential, int_space=5):
        self.basis_set = basis_set
        self.exchange_potential = exchange_potential
        self.correlation_potential = correlation_potential
        self.int_space = int_space
        self.density_matrix = np.matrix([])
        self.electron_density = {}
        self.integrand = {}

    def integrate(self, density_matrix, i, j):

        if density_matrix is not self.density_matrix:
            self.density_matrix = density_matrix
            self.electron_density = {}
            self.integrand = {}

        g_i = self.basis_set[i]
        g_j = self.basis_set[j]

        def electron_density(x, y, z):
            if (x, y, z) not in self.electron_density:
                density = 0
                for a, basis_a in enumerate(self.basis_set):
                    for b, basis_b in enumerate(self.basis_set):
                        if a == b:
                            density += density_matrix.item(a, b) * basis_a.value(x, y, z)**2
                        elif a <= b:
                            density += 2 * basis_a.value(x, y, z) * density_matrix.item(a, b) * basis_b.value(x, y, z)
                self.electron_density[(x, y, z)] = density
            return self.electron_density[(x, y, z)]

        def integrand(x, y, z):
            if (x, y, z) not in self.integrand:
                if i == j:
                    self.integrand[(x, y, z)] = g_i.value(x, y, z)**2 \
                    * (self.exchange_potential.calculate(electron_density(x, y, z))
                    + self.correlation_potential.calculate(electron_density(x, y, z)))
                else:
                    self.integrand[(x, y, z)] = g_i.value(x, y, z) * g_j.value(x, y, z) \
                    * (self.exchange_potential.calculate(electron_density(x, y, z))
                    + self.correlation_potential.calculate(electron_density(x, y, z)))
            return self.integrand[(x, y, z)]

        integral, error = integrate.nquad(integrand, [[-self.int_space, self.int_space],
        [-self.int_space, self.int_space], [-self.int_space, self.int_space]],
        opts=[{'epsabs': 1e-3, 'epsrel': 0}, {'epsabs': 1e-3, 'epsrel': 0}, {'epsabs': 1e-3, 'epsrel': 0}])

        return integral
