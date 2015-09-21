from unittest import TestCase
from unittest.mock import MagicMock
from unittest.mock import Mock
from src.main import Coulomb


class TestCoulombsLaw(TestCase):

    def setUp(self):
        self.num_nuc = 2
        self.mock_array = MagicMock()
        self.mock_nuclei = MagicMock()
        self.mock_array.__len__.return_value = self.num_nuc
        self.mock_array.__getitem__.return_value = self.mock_nuclei
        self.mock_nuclei.get_name.return_value = 'mockNuclei'
        self.coulomb = Coulomb(self.mock_array)

    def test_calculate_total_electric_potential_energy_callsLenMethod(self):
        self.coulomb.calculate_total_electric_potential_energy()
        self.mock_array.__len__.assert_called_once_with()

    def test_calculate_total_electric_potential_energy_callsGetNameMethod(self):
        self.coulomb.calculate_total_electric_potential_energy()
        self.assertEquals(self.mock_nuclei.get_name.call_count, (self.num_nuc * (self.num_nuc - 1)))

    def test_calculate_total_electric_potential_energy_returnsAFloat(self):
        energy = self.coulomb.calculate_total_electric_potential_energy()
        self.assertIsInstance(energy, float)
