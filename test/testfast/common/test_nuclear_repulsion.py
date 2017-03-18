from unittest import TestCase
from numpy import testing
from src.common import coulomb_matrix
from src.objects import Nuclei


class TestCoulomb(TestCase):

    def test_coulomb_calc_electric_potential_energy_returns_correct_potential_energy_of_hydrogen(self):
        nuclei_a = Nuclei('HYDROGEN', 1, 1, (0, 0, 0.70))
        nuclei_b = Nuclei('HYDROGEN', 1, 1, (0, 0, -0.7))
        nuclei_array = [nuclei_a, nuclei_b]
        energy = coulomb_matrix(nuclei_array)
        testing.assert_approx_equal(energy.item(0, 1), 0.7143, 4)

    def test_coulomb_calc_electric_potential_energy_returns_correct_potential_energy_of_heh_anion(self):
        nuclei_a = Nuclei('HYDROGEN', 1, 1, (0, 0, 0.7316))
        nuclei_b = Nuclei('HELIUM', 2, 4, (0, 0, -0.7316))
        nuclei_array = [nuclei_a, nuclei_b]
        energy = coulomb_matrix(nuclei_array)
        testing.assert_approx_equal(energy.item(0, 1), 1.366867, 7)
