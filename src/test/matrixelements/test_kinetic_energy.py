from unittest import TestCase
from src.main.matrixelements import KineticEnergyIntegral
from numpy import testing
from unittest.mock import MagicMock


class TestKineticEnergyIntegralHehCation(TestCase):

    def setUp(self):
        helium = MagicMock()
        hydrogen = MagicMock()
        basis_set_array = [helium, hydrogen]
        self.kinetic_energy_integral = KineticEnergyIntegral(basis_set_array)

        hydrogen.get_x.return_value = float(0)
        hydrogen.get_y.return_value = float(0)
        hydrogen.get_z.return_value = float(-0.7316)
        hydrogen.get_orbital_type.return_value = 'S'

        helium.get_x.return_value = float(0)
        helium.get_y.return_value = float(0)
        helium.get_z.return_value = float(0.7316)
        helium.get_orbital_type.return_value = 'S'

        hydrogen.get_array_of_coefficients.return_value = [[0.15432897, 3.42525091], [0.53532814, 0.62391373], [0.44463454, 0.16885540]]
        helium.get_array_of_coefficients.return_value = [[0.15432897, 9.75393461], [0.53532814, 1.77669115], [0.44463454, 0.48084429]]

    def test_method_calculate_returns_the_kinetic_energy_of_helium_1s(self):
        kinetic_energy = self.kinetic_energy_integral.calculate(0, 0)
        testing.assert_approx_equal(kinetic_energy, 2.1643, 5)

    def test_method_calculate_returs_the_kinetic_energy_of_quantum_mechanical_off_diagonal_parts(self):
        kinetic_energy = self.kinetic_energy_integral.calculate(1, 0)
        testing.assert_approx_equal(kinetic_energy, 0.167, 3)

    def test_method_calculate_returs_the_kinetic_energy_of_hydrogen_1s(self):
        kinetic_energy = self.kinetic_energy_integral.calculate(1, 1)
        testing.assert_approx_equal(kinetic_energy, 0.7600, 4)


class TestKineticEnergyIntegralHydrogen(TestCase):

    def setUp(self):
        hydrogen_1 = MagicMock()
        hydrogen_2 = MagicMock()
        basis_set_array = [hydrogen_1, hydrogen_2]
        self.kinetic_energy_integral = KineticEnergyIntegral(basis_set_array)

        hydrogen_1.get_x.return_value = float(0)
        hydrogen_1.get_y.return_value = float(0)
        hydrogen_1.get_z.return_value = float(0.7)
        hydrogen_1.get_orbital_type.return_value = 'S'

        hydrogen_2.get_x.return_value = float(0)
        hydrogen_2.get_y.return_value = float(0)
        hydrogen_2.get_z.return_value = float(-0.7)
        hydrogen_2.get_orbital_type.return_value = 'S'

        hydrogen_1.get_array_of_coefficients.return_value = [[0.15432897, 3.42525091], [0.53532814, 0.62391373], [0.44463454, 0.16885540]]
        hydrogen_2.get_array_of_coefficients.return_value = [[0.15432897, 3.42525091], [0.53532814, 0.62391373], [0.44463454, 0.16885540]]

    def test_method_calculate_returns_the_kinetic_energy_of_hydrogen_1s(self):
        kinetic_energy = self.kinetic_energy_integral.calculate(0, 0)
        testing.assert_approx_equal(kinetic_energy, 0.7600, 4)

    def test_method_calculate_returs_the_kinetic_energy_of_quantum_mechanical_off_diagonal_parts(self):
        kinetic_energy = self.kinetic_energy_integral.calculate(1, 0)
        testing.assert_approx_equal(kinetic_energy, 0.2365, 4)
