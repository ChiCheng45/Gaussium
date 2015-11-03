from unittest import TestCase
from numpy import testing
from unittest.mock import MagicMock
from src.main.matrixelements import KineticEnergyElement
import numpy as np


class TestKineticEnergyIntegralHehCation(TestCase):

    def setUp(self):
        helium = MagicMock()
        hydrogen = MagicMock()
        hydrogen_basis_1 = MagicMock()
        hydrogen_basis_2 = MagicMock()
        hydrogen_basis_3 = MagicMock()
        helium_basis_1 = MagicMock()
        helium_basis_2 = MagicMock()
        helium_basis_3 = MagicMock()
        basis_set_array = [helium, hydrogen]
        self.kinetic_energy_integral = KineticEnergyElement(basis_set_array)

        hydrogen.get_coordinates.return_value = np.matrix([[0], [0], [-0.7316]])
        helium.get_coordinates.return_value = np.matrix([[0], [0], [0.7316]])

        hydrogen.get_primitive_gaussian_array.return_value = [hydrogen_basis_1, hydrogen_basis_2, hydrogen_basis_3]
        helium.get_primitive_gaussian_array.return_value = [helium_basis_1, helium_basis_2, helium_basis_3]

        hydrogen_basis_1.get_contraction.return_value = 0.15432897
        hydrogen_basis_1.get_exponent.return_value = 3.42525091
        hydrogen_basis_1.get_orbital_type.return_value = 'S'
        hydrogen_basis_2.get_contraction.return_value = 0.53532814
        hydrogen_basis_2.get_exponent.return_value = 0.62391373
        hydrogen_basis_2.get_orbital_type.return_value = 'S'
        hydrogen_basis_3.get_contraction.return_value = 0.44463454
        hydrogen_basis_3.get_exponent.return_value = 0.16885540
        hydrogen_basis_3.get_orbital_type.return_value = 'S'

        helium_basis_1.get_contraction.return_value = 0.15432897
        helium_basis_1.get_exponent.return_value = 9.75393461
        helium_basis_1.get_orbital_type.return_value = 'S'
        helium_basis_2.get_contraction.return_value = 0.53532814
        helium_basis_2.get_exponent.return_value = 1.77669115
        helium_basis_2.get_orbital_type.return_value = 'S'
        helium_basis_3.get_contraction.return_value = 0.44463454
        helium_basis_3.get_exponent.return_value = 0.48084429
        helium_basis_3.get_orbital_type.return_value = 'S'

    def test_method_calculate_returns_the_kinetic_energy_of_helium_1s(self):
        kinetic_energy = self.kinetic_energy_integral.calculate(0, 0)
        testing.assert_approx_equal(kinetic_energy, 2.1643, 5)

    def test_method_calculate_returns_the_kinetic_energy_of_quantum_mechanical_off_diagonal_parts(self):
        kinetic_energy = self.kinetic_energy_integral.calculate(0, 1)
        testing.assert_approx_equal(kinetic_energy, 0.167, 3)

    def test_method_calculate_returns_the_kinetic_energy_of_hydrogen_1s(self):
        kinetic_energy = self.kinetic_energy_integral.calculate(1, 1)
        testing.assert_approx_equal(kinetic_energy, 0.7600, 4)


class TestKineticEnergyIntegralHydrogen(TestCase):

    def setUp(self):
        hydrogen_1 = MagicMock()
        hydrogen_2 = MagicMock()
        hydrogen_basis_1 = MagicMock()
        hydrogen_basis_2 = MagicMock()
        hydrogen_basis_3 = MagicMock()
        basis_set_array = [hydrogen_1, hydrogen_2]
        self.kinetic_energy_integral = KineticEnergyElement(basis_set_array)

        hydrogen_1.get_coordinates.return_value = np.matrix([[0], [0], [-0.7000]])
        hydrogen_2.get_coordinates.return_value = np.matrix([[0], [0], [0.7000]])

        hydrogen_1.get_primitive_gaussian_array.return_value = [hydrogen_basis_1, hydrogen_basis_2, hydrogen_basis_3]
        hydrogen_2.get_primitive_gaussian_array.return_value = [hydrogen_basis_1, hydrogen_basis_2, hydrogen_basis_3]

        hydrogen_basis_1.get_contraction.return_value = 0.15432897
        hydrogen_basis_1.get_exponent.return_value = 3.42525091
        hydrogen_basis_1.get_orbital_type.return_value = 'S'
        hydrogen_basis_2.get_contraction.return_value = 0.53532814
        hydrogen_basis_2.get_exponent.return_value = 0.62391373
        hydrogen_basis_2.get_orbital_type.return_value = 'S'
        hydrogen_basis_3.get_contraction.return_value = 0.44463454
        hydrogen_basis_3.get_exponent.return_value = 0.16885540
        hydrogen_basis_3.get_orbital_type.return_value = 'S'

    def test_method_calculate_returns_the_kinetic_energy_of_hydrogen_1s(self):
        kinetic_energy = self.kinetic_energy_integral.calculate(0, 0)
        testing.assert_approx_equal(kinetic_energy, 0.7600, 4)

    def test_method_calculate_returns_the_kinetic_energy_of_quantum_mechanical_off_diagonal_parts(self):
        kinetic_energy = self.kinetic_energy_integral.calculate(1, 0)
        testing.assert_approx_equal(kinetic_energy, 0.2365, 4)
