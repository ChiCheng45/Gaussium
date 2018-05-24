from scipy import integrate
import numpy as np
import quadpy


class ExchangeCorrelation:

    def __init__(self, basis_set, exchange_potential, correlation_potential, int_space=25, epsabs=1e-9, epsrel=0.0):
        self.basis_set = basis_set
        self.exchange_potential = exchange_potential
        self.correlation_potential = correlation_potential
        self.int_space = int_space
        self.epsabs = epsabs
        self.epsrel = epsrel

    def integrate(self, density_matrix, i, j):

        def electron_density(x, y, z):
            density = 0
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
            return self.basis_set[i].value(x, y, z) * (self.exchange_potential.calculate(electron_density(x, y, z))
            + self.correlation_potential.calculate(electron_density(x, y, z))) * self.basis_set[j].value(x, y, z) \
            * rho**2

        def integrand_spherical(rho):
            return quadpy.sphere.integrate_spherical(
                lambda azimuthal, polar: integrand(rho, polar, azimuthal),
                rule=quadpy.sphere.Lebedev(3)
            )

        return integrate.quad(integrand_spherical, 0.0, self.int_space, epsabs=self.epsabs, epsrel=self.epsrel)[0]