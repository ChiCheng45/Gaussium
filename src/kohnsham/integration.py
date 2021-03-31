from scipy import integrate
import numpy as np
import quadpy


class ExchangeCorrelation:

    def __init__(self, basis_set, exchange, correlation, int_space=25, epsabs=1e-9, epsrel=0.0):
        self.basis_set = basis_set
        self.exchange = exchange
        self.correlation = correlation
        self.int_space = int_space
        self.epsabs = epsabs
        self.epsrel = epsrel

    def integrate_potential(self, density_matrix, i, j):

        g_1 = self.basis_set[i]
        g_2 = self.basis_set[j]

        def electron_density(x, y, z):
            density = 0.0
            for a, basis_a in enumerate(self.basis_set):
                for b, basis_b in enumerate(self.basis_set):
                    if a == b:
                        density += density_matrix.item(a, b) * basis_a.value(x, y, z)**2
                    elif a < b:
                        density += 2 * basis_a.value(x, y, z) * density_matrix.item(a, b) * basis_b.value(x, y, z)
            return density

        def integrand(rho, theta, phi):
            x = rho * np.sin(theta) * np.cos(phi)
            y = rho * np.sin(theta) * np.sin(phi)
            z = rho * np.cos(theta)
            density = electron_density(x, y, z)
            return g_1.value(x, y, z) * (self.exchange.potential(density)
            + self.correlation.potential(density)) * g_2.value(x, y, z) \
            * rho**2

        def integrand_spherical(rho):
            scheme = quadpy.u3.schemes["lebedev_065"]()
            return scheme.integrate_spherical(
                lambda angles: integrand(rho, angles[0], angles[1])
            )

        return integrate.quad(
            integrand_spherical, 0.0, self.int_space, epsabs=self.epsabs,
            epsrel=self.epsrel, points=[0.0]
        )[0]

    def integrate_energy(self, density_matrix):

        def electron_density(x, y, z):
            density = 0.0
            for a, basis_a in enumerate(self.basis_set):
                for b, basis_b in enumerate(self.basis_set):
                    if a == b:
                        density += density_matrix.item(a, b) * basis_a.value(x, y, z)**2
                    elif a < b:
                        density += 2 * basis_a.value(x, y, z) * density_matrix.item(a, b) * basis_b.value(x, y, z)
            return density

        def integrand(rho, theta, phi):
            x = rho * np.sin(theta) * np.cos(phi)
            y = rho * np.sin(theta) * np.sin(phi)
            z = rho * np.cos(theta)
            density = electron_density(x, y, z)
            return density * (self.exchange.energy(density)
            + self.correlation.energy(density)) * rho**2

        def integrand_spherical(rho):
            scheme = quadpy.u3.schemes["lebedev_065"]()
            return scheme.integrate_spherical(
                lambda angles: integrand(rho, angles[0], angles[1])
            )

        return integrate.quad(
            integrand_spherical, 0.0, self.int_space, epsabs=self.epsabs,
            epsrel=self.epsrel, points=[0.0]
        )[0]
