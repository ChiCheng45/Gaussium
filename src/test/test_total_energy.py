from unittest import TestCase
import numpy as np
from numpy import testing
from src.main.self_consistent_field import TotalEnergy


class TestTotalEnergyHeH(TestCase):

    def setUp(self):
        self.density_matrix = np.matrix([[1.7266, 0.2599], [0.2599, 0.0391]])
        self.h_core_matrix = np.matrix([[-2.6527, -1.3472], [-1.3472, -1.7318]])
        self.fock_matrix = np.matrix([[-1.3904, -0.9732], [-0.9732, -0.7429]])

    def test_calculate_total_energy_for_the_heh_case(self):
        total_energy = TotalEnergy()
        electronic_energy = total_energy.calculate_total_energy(self.density_matrix, self.h_core_matrix, self.fock_matrix)
        testing.assert_approx_equal(electronic_energy, -4.14186, 6)
