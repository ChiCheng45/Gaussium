from unittest import TestCase
from unittest.mock import MagicMock
from src.main import TwoElectronRepulsion


class TestTwoElectronRepulsion(TestCase):

    def setUp(self):
        self.mock_matrix = MagicMock()
        self.mock_array = MagicMock()
        self.mock_filereader = MagicMock()
        self.mock_int = MagicMock()
        self.mock_nuclei = MagicMock()
        self.mock_array.__getitem__.return_value = self.mock_nuclei
        self.twoelectronrepulsion = TwoElectronRepulsion()

    def test_calculate_integral_returns_null(self):
        energy = self.twoelectronrepulsion.calculate_integral(self.mock_array, self.mock_filereader, self.mock_int, self.mock_int)
        self.assertEquals(0, energy)

    def test_calculate_integral_calls_create_basis_set_array(self):
        self.twoelectronrepulsion.calculate_integral(self.mock_array, self.mock_filereader, self.mock_int, self.mock_int)
        self.assertEquals(self.mock_filereader.create_basis_set_array.call_count, 2)

    def test_calculate_integral_calls_get_name_method_of_mock_nuclei(self):
        self.twoelectronrepulsion.calculate_integral(self.mock_array, self.mock_filereader, self.mock_int, self.mock_int)
        self.assertEquals(self.mock_nuclei.get_name.call_count, 2)

