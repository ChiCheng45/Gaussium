from unittest import TestCase
from src.main.common import CoulombsLaw
from src.main.objects import Nuclei
from numpy import testing


class TestCoulomb(TestCase):

    def test_coulomb_calc_electric_potential_energy_returns_correct_potential_energy_of_hydrogen(self):
        nuclei_a = Nuclei(['HYDROGEN', 1, 1, 0, 0, 0.7])
        nuclei_b = Nuclei(['HYDROGEN', 1, 1, 0, 0, -0.7])
        energy = CoulombsLaw(nuclei_a, nuclei_b).calc_electric_potential_energy()
        testing.assert_approx_equal(energy, 0.7143, 4)

    def test_coulomb_calc_electric_potential_energy_returns_correct_potential_energy_of_heh_anion(self):
        nuclei_a = Nuclei(['HYDROGEN', 1, 1, 0, 0, 0.7316])
        nuclei_b = Nuclei(['HELIUM', 2, 4, 0, 0, -0.7316])
        energy = CoulombsLaw(nuclei_a, nuclei_b).calc_electric_potential_energy()
        testing.assert_approx_equal(energy, 1.366867, 7)
