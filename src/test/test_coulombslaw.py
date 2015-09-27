from unittest import TestCase
from unittest.mock import MagicMock
from src.main import CoulombTotal
from src.main import Coulomb
from src.main import Nuclei
from numpy import testing


class TestCoulombTotal(TestCase):

    def setUp(self):
        self.num = 11
        self.mock_array = MagicMock()
        self.mock_nuclei = MagicMock()
        self.mock_coulomb = MagicMock()
        self.mock_energy = MagicMock()
        self.mock_array.__len__.return_value = self.num
        self.mock_array.__getitem__.return_value = self.mock_nuclei
        self.coulomb_total = CoulombTotal(self.mock_coulomb, self.mock_array)

    def test_coulomb_total_calculate_total_electric_potential_energy_calls_len_method_on_mock_array(self):
        self.coulomb_total.calculate_total_electric_potential_energy()
        self.mock_array.__len__.assert_called_once_with()

    def test_coulomb_total_calculate_total_electric_potential_energy_calls_get_item_method_on_mock_array(self):
        self.coulomb_total.calculate_total_electric_potential_energy()
        self.assertEquals(self.mock_array.__getitem__.call_count, 220)

    def test_coulomb_total_calculate_total_electric_potential_energy_calls_calc_electric_potential_energy_on_mock_coulomb_object(self):
        self.coulomb_total.calculate_total_electric_potential_energy()
        self.assertEquals(
            self.mock_coulomb(self.mock_nuclei, self.mock_nuclei).calc_electric_potential_energy.call_count,
            110)

    def test_coulomb_total_calculate_total_electric_potential_energy_returns_a_n_by_n_matrix(self):
        matrix = self.coulomb_total.calculate_total_electric_potential_energy()
        self.assertEquals(matrix.shape, (11, 11))

    def test_coulomb_total_calculate_total_electric_potential_energy_returns_a_symmetric_matrix(self):
        matrix = self.coulomb_total.calculate_total_electric_potential_energy()
        testing.assert_array_equal(matrix, matrix.T)

    def test_coulomb_total_calculate_total_electric_potential_energy_returns_a_matrix_with_trace_eq_0(self):
        matrix = self.coulomb_total.calculate_total_electric_potential_energy()
        self.assertEquals(matrix.trace(), 0)

    def test_coulomb_total_calculate_total_electric_potential_energy_correct_potential_energy_matrix_of_hydrogen(self):
        nuclei_a = Nuclei(['HYDROGEN', 1, 1, 0, 0, 0.7])
        nuclei_b = Nuclei(['HYDROGEN', 1, 1, 0, 0, -0.7])
        nuclei_array = [nuclei_a, nuclei_b]
        matrix = CoulombTotal(Coulomb, nuclei_array).calculate_total_electric_potential_energy()
        testing.assert_approx_equal(matrix.sum() / 2, 0.7143, 4)

    def test_coulomb_calc_electric_potential_energy_returns_correct_potential_energy_of_heh_anion(self):
        nuclei_a = Nuclei(['HYDROGEN', 1, 1, 0, 0, 0.7316])
        nuclei_b = Nuclei(['HELIUM', 2, 4, 0, 0, -0.7316])
        nuclei_array = [nuclei_a, nuclei_b]
        matrix = CoulombTotal(Coulomb, nuclei_array).calculate_total_electric_potential_energy()
        testing.assert_approx_equal(matrix.sum() / 2, 1.366867, 7)


class TestCoulomb(TestCase):

    def setUp(self):
        self.mock_nuclei = MagicMock()
        self.coulomb = Coulomb(self.mock_nuclei, self.mock_nuclei)

    def test_coulomb_calc_electric_potential_energy_calls_get_x_methods(self):
        self.coulomb.calc_electric_potential_energy()
        self.assertEquals(self.mock_nuclei.get_x.call_count, 2)

    def test_coulomb_calc_electric_potential_energy_calls_get_y_methods(self):
        self.coulomb.calc_electric_potential_energy()
        self.assertEquals(self.mock_nuclei.get_y.call_count, 2)

    def test_coulomb_calc_electric_potential_energy_calls_get_z_methods(self):
        self.coulomb.calc_electric_potential_energy()
        self.assertEquals(self.mock_nuclei.get_z.call_count, 2)

    def test_coulomb_calc_electric_potential_energy_returns_correct_potential_energy_of_hydrogen(self):
        nuclei_a = Nuclei(['HYDROGEN', 1, 1, 0, 0, 0.7])
        nuclei_b = Nuclei(['HYDROGEN', 1, 1, 0, 0, -0.7])
        energy = Coulomb(nuclei_a, nuclei_b).calc_electric_potential_energy()
        testing.assert_approx_equal(energy, 0.7143, 4)

    def test_coulomb_calc_electric_potential_energy_returns_correct_potential_energy_of_heh_anion(self):
        nuclei_a = Nuclei(['HYDROGEN', 1, 1, 0, 0, 0.7316])
        nuclei_b = Nuclei(['HELIUM', 2, 4, 0, 0, -0.7316])
        energy = Coulomb(nuclei_a, nuclei_b).calc_electric_potential_energy()
        testing.assert_approx_equal(energy, 1.366867, 7)
