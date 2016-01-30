from unittest import TestCase
from unittest.mock import MagicMock
from src.main.matrixelements import OrbitalOverlapMatrix
from numpy import testing


class TestOrbitalOverlapMatrixHeHCation(TestCase):

    def setUp(self):
        hydrogen_basis_1 = MagicMock(contraction=0.15432897, exponent=3.42525091, coordinates=(0, 0, -0.7316), integral_exponents=(0, 0, 0), normalisation = 1.794441832218435)
        hydrogen_basis_2 = MagicMock(contraction=0.53532814, exponent=0.62391373, coordinates=(0, 0, -0.7316), integral_exponents=(0, 0, 0), normalisation = 0.5003264923314032)
        hydrogen_basis_3 = MagicMock(contraction=0.44463454, exponent=0.16885540, coordinates=(0, 0, -0.7316), integral_exponents=(0, 0, 0), normalisation = 0.18773545851092535)
        helium_basis_1 = MagicMock(contraction=0.15432897, exponent=9.75393461, coordinates=(0, 0, 0.7316), integral_exponents=(0, 0, 0), normalisation = 3.9336432656254527)
        helium_basis_2 = MagicMock(contraction=0.53532814, exponent=1.77669115, coordinates=(0, 0, 0.7316), integral_exponents=(0, 0, 0), normalisation = 1.0967787981767012)
        helium_basis_3 = MagicMock(contraction=0.44463454, exponent=0.48084429, coordinates=(0, 0, 0.7316), integral_exponents=(0, 0, 0), normalisation = 0.41154131374122654)
        helium = MagicMock(primitive_gaussian_array=[helium_basis_1, helium_basis_2, helium_basis_3], coordinates=(0, 0, -0.7316), integral_exponents=(0, 0, 0))
        hydrogen = MagicMock(primitive_gaussian_array=[hydrogen_basis_1, hydrogen_basis_2, hydrogen_basis_3], coordinates=(0, 0, 0.7316), integral_exponents=(0, 0, 0))
        basis_set_array = [helium, hydrogen]

        self.kinetic_energy_integral = OrbitalOverlapMatrix(basis_set_array)

    def test_method_calculate_returns_the_overlap_of_he_1s_and_h_1s(self):
        kinetic_energy = self.kinetic_energy_integral.calculate(0, 1)
        testing.assert_approx_equal(kinetic_energy, 0.4508, 3)


class TestOrbitalOverlapMatrixHydrogen(TestCase):

    def setUp(self):
        hydrogen_basis_1 = MagicMock(contraction=0.15432897, exponent=3.42525091, coordinates=(0, 0, 0.7000), integral_exponents=(0, 0, 0), normalisation = 1.794441832218435)
        hydrogen_basis_2 = MagicMock(contraction=0.53532814, exponent=0.62391373, coordinates=(0, 0, 0.7000), integral_exponents=(0, 0, 0), normalisation = 0.5003264923314032)
        hydrogen_basis_3 = MagicMock(contraction=0.44463454, exponent=0.16885540, coordinates=(0, 0, 0.7000), integral_exponents=(0, 0, 0), normalisation = 0.18773545851092535)
        hydrogen_basis_4 = MagicMock(contraction=0.15432897, exponent=3.42525091, coordinates=(0, 0, -0.7000), integral_exponents=(0, 0, 0), normalisation = 1.794441832218435)
        hydrogen_basis_5 = MagicMock(contraction=0.53532814, exponent=0.62391373, coordinates=(0, 0, -0.7000), integral_exponents=(0, 0, 0), normalisation = 0.5003264923314032)
        hydrogen_basis_6 = MagicMock(contraction=0.44463454, exponent=0.16885540, coordinates=(0, 0, -0.7000), integral_exponents=(0, 0, 0), normalisation = 0.18773545851092535)
        hydrogen_1 = MagicMock(primitive_gaussian_array=[hydrogen_basis_1, hydrogen_basis_2, hydrogen_basis_3], coordinates=(0, 0, -0.7000), integral_exponents=(0, 0, 0))
        hydrogen_2 = MagicMock(primitive_gaussian_array=[hydrogen_basis_4, hydrogen_basis_5, hydrogen_basis_6], coordinates=(0, 0, 0.7000), integral_exponents=(0, 0, 0))
        basis_set_array = [hydrogen_1, hydrogen_2]

        self.kinetic_energy_integral = OrbitalOverlapMatrix(basis_set_array)

    def test_method_calculate_returns_the_overlap_of_h1_1s_and_h2_1s(self):
        kinetic_energy = self.kinetic_energy_integral.calculate(0, 1)
        testing.assert_approx_equal(kinetic_energy, 0.6593, 4)
