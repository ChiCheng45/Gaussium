from unittest import TestCase
from unittest.mock import MagicMock
from src.main.matrixelements import TwoElectronRepulsionElementCook, TwoElectronRepulsionElementOS, TwoElectronRepulsionElementHGP
from numpy import testing


class TestTwoElectronRepulsionElementCook(TestCase):

    def setUp(self):
        hydrogen_basis_1 = MagicMock(contraction=0.15432897, exponent=3.42525091, coordinates=(0, 0, -0.7316), integral_exponents=(0, 0, 0))
        hydrogen_basis_2 = MagicMock(contraction=0.53532814, exponent=0.62391373, coordinates=(0, 0, -0.7316), integral_exponents=(0, 0, 0))
        hydrogen_basis_3 = MagicMock(contraction=0.44463454, exponent=0.16885540, coordinates=(0, 0, -0.7316), integral_exponents=(0, 0, 0))
        helium_basis_1 = MagicMock(contraction=0.15432897, exponent=9.75393461, coordinates=(0, 0, 0.7316), integral_exponents=(0, 0, 0))
        helium_basis_2 = MagicMock(contraction=0.53532814, exponent=1.77669115, coordinates=(0, 0, 0.7316), integral_exponents=(0, 0, 0))
        helium_basis_3 = MagicMock(contraction=0.44463454, exponent=0.48084429, coordinates=(0, 0, 0.7316), integral_exponents=(0, 0, 0))
        hydrogen_basis_1.normalisation.return_value = 1.794441832218435
        hydrogen_basis_2.normalisation.return_value = 0.5003264923314032
        hydrogen_basis_3.normalisation.return_value = 0.18773545851092535
        helium_basis_1.normalisation.return_value = 3.9336432656254527
        helium_basis_2.normalisation.return_value = 1.0967787981767012
        helium_basis_3.normalisation.return_value = 0.41154131374122654
        helium = MagicMock(primitive_gaussian_array=[helium_basis_1, helium_basis_2, helium_basis_3], coordinates=(0, 0, 0.7316), integral_exponents=(0, 0, 0))
        hydrogen = MagicMock(primitive_gaussian_array=[hydrogen_basis_1, hydrogen_basis_2, hydrogen_basis_3], coordinates=(0, 0, -0.7316), integral_exponents=(0, 0, 0))
        basis_set_array = [helium, hydrogen]
        self.two_electron_repulsion = TwoElectronRepulsionElementCook(basis_set_array)

    def test_method_calculate_returns_element_for_0000(self):
        element = self.two_electron_repulsion.calculate(0, 0, 0, 0)
        testing.assert_approx_equal(element, 1.3072, 4)

    def test_method_calculate_returns_element_for_1000(self):
        element = self.two_electron_repulsion.calculate(1, 0, 0, 0)
        testing.assert_approx_equal(element, 0.4373, 4)

    def test_method_calculate_returns_element_for_1010(self):
        element = self.two_electron_repulsion.calculate(1, 0, 1, 0)
        testing.assert_approx_equal(element, 0.1773, 4)

    def test_method_calculate_returns_element_for_1100(self):
        element = self.two_electron_repulsion.calculate(1, 1, 0, 0)
        testing.assert_approx_equal(element, 0.6057, 4)

    def test_method_calculate_returns_element_for_1110(self):
        element = self.two_electron_repulsion.calculate(1, 1, 1, 0)
        testing.assert_approx_equal(element, 0.3118, 4)

    def test_method_calculate_returns_element_for_1111(self):
        element = self.two_electron_repulsion.calculate(1, 1, 1, 1)
        testing.assert_approx_equal(element, 0.7746, 4)


class TestTwoElectronRepulsionElementOS(TestCase):

    def setUp(self):
        hydrogen_basis_1 = MagicMock(contraction=0.15432897, exponent=3.42525091, coordinates=(0, 0, -0.7316), integral_exponents=(0, 0, 0))
        hydrogen_basis_2 = MagicMock(contraction=0.53532814, exponent=0.62391373, coordinates=(0, 0, -0.7316), integral_exponents=(0, 0, 0))
        hydrogen_basis_3 = MagicMock(contraction=0.44463454, exponent=0.16885540, coordinates=(0, 0, -0.7316), integral_exponents=(0, 0, 0))
        helium_basis_1 = MagicMock(contraction=0.15432897, exponent=9.75393461, coordinates=(0, 0, 0.7316), integral_exponents=(0, 0, 0))
        helium_basis_2 = MagicMock(contraction=0.53532814, exponent=1.77669115, coordinates=(0, 0, 0.7316), integral_exponents=(0, 0, 0))
        helium_basis_3 = MagicMock(contraction=0.44463454, exponent=0.48084429, coordinates=(0, 0, 0.7316), integral_exponents=(0, 0, 0))
        hydrogen_basis_1.normalisation.return_value = 1.794441832218435
        hydrogen_basis_2.normalisation.return_value = 0.5003264923314032
        hydrogen_basis_3.normalisation.return_value = 0.18773545851092535
        helium_basis_1.normalisation.return_value = 3.9336432656254527
        helium_basis_2.normalisation.return_value = 1.0967787981767012
        helium_basis_3.normalisation.return_value = 0.41154131374122654
        helium = MagicMock(primitive_gaussian_array=[helium_basis_1, helium_basis_2, helium_basis_3], coordinates=(0, 0, 0.7316), integral_exponents=(0, 0, 0))
        hydrogen = MagicMock(primitive_gaussian_array=[hydrogen_basis_1, hydrogen_basis_2, hydrogen_basis_3], coordinates=(0, 0, -0.7316), integral_exponents=(0, 0, 0))
        basis_set_array = [helium, hydrogen]
        self.two_electron_repulsion = TwoElectronRepulsionElementOS(basis_set_array)

    def test_method_calculate_returns_element_for_0000(self):
        element = self.two_electron_repulsion.calculate(0, 0, 0, 0)
        testing.assert_approx_equal(element, 1.3072, 4)

    def test_method_calculate_returns_element_for_1000(self):
        element = self.two_electron_repulsion.calculate(1, 0, 0, 0)
        testing.assert_approx_equal(element, 0.4373, 4)

    def test_method_calculate_returns_element_for_1010(self):
        element = self.two_electron_repulsion.calculate(1, 0, 1, 0)
        testing.assert_approx_equal(element, 0.1773, 4)

    def test_method_calculate_returns_element_for_1100(self):
        element = self.two_electron_repulsion.calculate(1, 1, 0, 0)
        testing.assert_approx_equal(element, 0.6057, 4)

    def test_method_calculate_returns_element_for_1110(self):
        element = self.two_electron_repulsion.calculate(1, 1, 1, 0)
        testing.assert_approx_equal(element, 0.3118, 4)

    def test_method_calculate_returns_element_for_1111(self):
        element = self.two_electron_repulsion.calculate(1, 1, 1, 1)
        testing.assert_approx_equal(element, 0.7746, 4)


class TestTwoElectronRepulsionElementHGP(TestCase):

    def setUp(self):
        hydrogen_basis_1 = MagicMock(contraction=0.15432897, exponent=3.42525091, coordinates=(0, 0, -0.7316), integral_exponents=(0, 0, 0))
        hydrogen_basis_2 = MagicMock(contraction=0.53532814, exponent=0.62391373, coordinates=(0, 0, -0.7316), integral_exponents=(0, 0, 0))
        hydrogen_basis_3 = MagicMock(contraction=0.44463454, exponent=0.16885540, coordinates=(0, 0, -0.7316), integral_exponents=(0, 0, 0))
        helium_basis_1 = MagicMock(contraction=0.15432897, exponent=9.75393461, coordinates=(0, 0, 0.7316), integral_exponents=(0, 0, 0))
        helium_basis_2 = MagicMock(contraction=0.53532814, exponent=1.77669115, coordinates=(0, 0, 0.7316), integral_exponents=(0, 0, 0))
        helium_basis_3 = MagicMock(contraction=0.44463454, exponent=0.48084429, coordinates=(0, 0, 0.7316), integral_exponents=(0, 0, 0))
        hydrogen_basis_1.normalisation.return_value = 1.794441832218435
        hydrogen_basis_2.normalisation.return_value = 0.5003264923314032
        hydrogen_basis_3.normalisation.return_value = 0.18773545851092535
        helium_basis_1.normalisation.return_value = 3.9336432656254527
        helium_basis_2.normalisation.return_value = 1.0967787981767012
        helium_basis_3.normalisation.return_value = 0.41154131374122654
        helium = MagicMock(primitive_gaussian_array=[helium_basis_1, helium_basis_2, helium_basis_3], coordinates=(0, 0, 0.7316), integral_exponents=(0, 0, 0))
        hydrogen = MagicMock(primitive_gaussian_array=[hydrogen_basis_1, hydrogen_basis_2, hydrogen_basis_3], coordinates=(0, 0, -0.7316), integral_exponents=(0, 0, 0))
        basis_set_array = [helium, hydrogen]
        self.two_electron_repulsion = TwoElectronRepulsionElementHGP(basis_set_array)

    def test_method_calculate_returns_element_for_0000(self):
        element = self.two_electron_repulsion.calculate(0, 0, 0, 0)
        testing.assert_approx_equal(element, 1.3072, 4)

    def test_method_calculate_returns_element_for_1000(self):
        element = self.two_electron_repulsion.calculate(1, 0, 0, 0)
        testing.assert_approx_equal(element, 0.4373, 4)

    def test_method_calculate_returns_element_for_1010(self):
        element = self.two_electron_repulsion.calculate(1, 0, 1, 0)
        testing.assert_approx_equal(element, 0.1773, 4)

    def test_method_calculate_returns_element_for_1100(self):
        element = self.two_electron_repulsion.calculate(1, 1, 0, 0)
        testing.assert_approx_equal(element, 0.6057, 4)

    def test_method_calculate_returns_element_for_1110(self):
        element = self.two_electron_repulsion.calculate(1, 1, 1, 0)
        testing.assert_approx_equal(element, 0.3118, 4)

    def test_method_calculate_returns_element_for_1111(self):
        element = self.two_electron_repulsion.calculate(1, 1, 1, 1)
        testing.assert_approx_equal(element, 0.7746, 4)
