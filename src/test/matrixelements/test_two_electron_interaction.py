from unittest import TestCase
from unittest.mock import MagicMock
from src.main.matrixelements import TwoElectronRepulsionElement
from numpy import testing


class TestTwoElectronRepulsionElement(TestCase):

    def setUp(self):
        hydrogen_basis_1 = MagicMock(contraction=0.15432897, exponent=3.42525091, coordinates=(0.000000, 0.000000, -0.7316), integral_exponents=(0, 0, 0))
        hydrogen_basis_2 = MagicMock(contraction=0.53532814, exponent=0.44463454, coordinates=(0.000000, 0.000000, -0.7316), integral_exponents=(0, 0, 0))
        hydrogen_basis_3 = MagicMock(contraction=0.44463454, exponent=0.16885540, coordinates=(0.000000, 0.000000, -0.7316), integral_exponents=(0, 0, 0))
        helium_basis_1 = MagicMock(contraction=0.15432897, exponent=9.75393461, coordinates=(0.000000, 0.000000, 0.7316), integral_exponents=(0, 0, 0))
        helium_basis_2 = MagicMock(contraction=0.53532814, exponent=1.77669115, coordinates=(0.000000, 0.000000, 0.7316), integral_exponents=(0, 0, 0))
        helium_basis_3 = MagicMock(contraction=0.44463454, exponent=0.48084429, coordinates=(0.000000, 0.000000, 0.7316), integral_exponents=(0, 0, 0))
        hydrogen = MagicMock(primitive_gaussian_array=[hydrogen_basis_1, hydrogen_basis_2, hydrogen_basis_3])
        helium = MagicMock(primitive_gaussian_array=[helium_basis_1, helium_basis_2, helium_basis_3])

        basis_set_array = [helium, hydrogen]
        self.two_electron_repulsion = TwoElectronRepulsionElement(basis_set_array)

    def test_method_calculate_returns_element_for_0000(self):
        element = self.two_electron_repulsion.calculate(0, 0, 0, 0)
        testing.assert_approx_equal(element, 1.3072, 4)

    def test_method_calculate_returns_element_for_0001(self):
        element = self.two_electron_repulsion.calculate(0, 1, 0, 0)
        print(element)
        testing.assert_approx_equal(element, 1.3072, 4)

    def test_store_integrals(self):
        dictionary = self.two_electron_repulsion.store_integrals()
        print(dictionary)
