import itertools


class KSEnergy:

    def __init__(self, number_electrons, exchange_correlation):
        self.number_electrons = number_electrons
        self.exchange_correlation = exchange_correlation

    def restricted(self, orbital_energies, density_matrix, elst_matrix, xc_matrix):
        total_energy = 2 * sum(orbital_energies[:self.number_electrons // 2])
        length = density_matrix.shape[0]
        for i, j in itertools.product(range(length), repeat=2):
            total_energy -= density_matrix.item(i, j) * (
                    (1/2) * elst_matrix.item(i, j) + xc_matrix.item(i, j)
            )
        total_energy += self.exchange_correlation.integrate_energy(density_matrix)
        return total_energy
