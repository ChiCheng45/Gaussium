from unittest import TestCase
import numpy as np
from numpy import testing
from gaussium.hartreefock import TotalEnergy


class TestTotalEnergyHeH(TestCase):

    def setUp(self):
        self.density_matrix = np.array([
                [1.7266, 0.2599],
                [0.2599, 0.0391]
        ])

        self.h_core_matrix = np.array([
                [-2.6527, -1.3472],
                [-1.3472, -1.7318]
        ])

        self.fock_matrix = np.array([
                [-1.3904, -0.9732],
                [-0.9732, -0.7429]
        ])

    def test_calculate_total_energy_for_the_heh_case(self):
        total_energy = TotalEnergy(self.h_core_matrix)
        electronic_energy = total_energy.restricted(self.density_matrix, self.fock_matrix)
        testing.assert_approx_equal(electronic_energy, -4.14186, 6)
