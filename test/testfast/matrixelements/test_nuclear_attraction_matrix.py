from unittest import TestCase
from unittest.mock import MagicMock

from numpy import testing

from src.matrixelements import NuclearAttractionMatrix


class TestNuclearAttractionMatrixHeHCation(TestCase):

    def setUp(self):
        hydrogen_basis_1 = MagicMock(contraction=0.15432897, exponent=3.42525091, coordinates=(0, 0, -0.7316), integral_exponents=(0, 0, 0), normalisation = 1.794441832218435)
        hydrogen_basis_2 = MagicMock(contraction=0.53532814, exponent=0.62391373, coordinates=(0, 0, -0.7316), integral_exponents=(0, 0, 0), normalisation = 0.5003264923314032)
        hydrogen_basis_3 = MagicMock(contraction=0.44463454, exponent=0.16885540, coordinates=(0, 0, -0.7316), integral_exponents=(0, 0, 0), normalisation = 0.18773545851092535)
        helium_basis_1 = MagicMock(contraction=0.15432897, exponent=9.75393461, coordinates=(0, 0, 0.7316), integral_exponents=(0, 0, 0), normalisation = 3.9336432656254527)
        helium_basis_2 = MagicMock(contraction=0.53532814, exponent=1.77669115, coordinates=(0, 0, 0.7316), integral_exponents=(0, 0, 0), normalisation = 1.0967787981767012)
        helium_basis_3 = MagicMock(contraction=0.44463454, exponent=0.48084429, coordinates=(0, 0, 0.7316), integral_exponents=(0, 0, 0), normalisation = 0.41154131374122654)
        helium_basis = MagicMock(primitive_gaussian_array=[helium_basis_1, helium_basis_2, helium_basis_3], coordinates=(0, 0, -0.7316), integral_exponents=(0, 0, 0))
        hydrogen_basis = MagicMock(primitive_gaussian_array=[hydrogen_basis_1, hydrogen_basis_2, hydrogen_basis_3], coordinates=(0, 0, 0.7316), integral_exponents=(0, 0, 0))
        basis_set_array = [helium_basis, hydrogen_basis]

        helium_nuclei = MagicMock(element='HELIUM', charge=2, mass=4, coordinates=(0.000000, 0.000000, 0.7316))
        hydrogen_nuclei = MagicMock(element='HYDROGEN', charge=1, mass=1, coordinates=(0.000000, 0.000000, -0.7316))
        nuclei_array = [helium_nuclei, hydrogen_nuclei]

        self.kinetic_energy_integral = NuclearAttractionMatrix(basis_set_array, nuclei_array)

    def test_method_calculate_returns_the_nuclear_attraction_energy_of_helium_1s(self):
        kinetic_energy = self.kinetic_energy_integral.calculate(0, 0)
        testing.assert_approx_equal(kinetic_energy, -4.817, 5)

    def test_method_calculate_returns_the_nuclear_attraction_energy_of_quantum_mechanical_off_diagonal_parts(self):
        kinetic_energy = self.kinetic_energy_integral.calculate(0, 1)
        testing.assert_approx_equal(kinetic_energy, -1.5142, 3)

    def test_method_calculate_returns_the_nuclear_attraction_energy_of_hydrogen_1s(self):
        kinetic_energy = self.kinetic_energy_integral.calculate(1, 1)
        testing.assert_approx_equal(kinetic_energy, -2.4912, 4)


class TestNuclearAttractionMatrixHydrogen(TestCase):

    def setUp(self):
        hydrogen_basis_1 = MagicMock(contraction=0.15432897, exponent=3.42525091, coordinates=(0, 0, 0.7000), integral_exponents=(0, 0, 0), normalisation = 1.794441832218435)
        hydrogen_basis_2 = MagicMock(contraction=0.53532814, exponent=0.62391373, coordinates=(0, 0, 0.7000), integral_exponents=(0, 0, 0), normalisation = 0.5003264923314032)
        hydrogen_basis_3 = MagicMock(contraction=0.44463454, exponent=0.16885540, coordinates=(0, 0, 0.7000), integral_exponents=(0, 0, 0), normalisation = 0.18773545851092535)
        hydrogen_basis_4 = MagicMock(contraction=0.15432897, exponent=3.42525091, coordinates=(0, 0, -0.7000), integral_exponents=(0, 0, 0), normalisation = 1.794441832218435)
        hydrogen_basis_5 = MagicMock(contraction=0.53532814, exponent=0.62391373, coordinates=(0, 0, -0.7000), integral_exponents=(0, 0, 0), normalisation = 0.5003264923314032)
        hydrogen_basis_6 = MagicMock(contraction=0.44463454, exponent=0.16885540, coordinates=(0, 0, -0.7000), integral_exponents=(0, 0, 0), normalisation = 0.18773545851092535)
        hydrogen_1_basis = MagicMock(primitive_gaussian_array=[hydrogen_basis_1, hydrogen_basis_2, hydrogen_basis_3], coordinates=(0, 0, -0.7000), integral_exponents=(0, 0, 0))
        hydrogen_2_basis = MagicMock(primitive_gaussian_array=[hydrogen_basis_4, hydrogen_basis_5, hydrogen_basis_6], coordinates=(0, 0, 0.7000), integral_exponents=(0, 0, 0))
        basis_set_array = [hydrogen_1_basis, hydrogen_2_basis]

        hydrogen_1_nuclei = MagicMock(element='HYDROGEN', charge=1, mass=1, coordinates=(0.000000, 0.000000, 0.7000))
        hydrogen_2_nuclei = MagicMock(element='HYDROGEN', charge=1, mass=1, coordinates=(0.000000, 0.000000, -0.7000))
        nuclei_array = [hydrogen_1_nuclei, hydrogen_2_nuclei]

        self.kinetic_energy_integral = NuclearAttractionMatrix(basis_set_array, nuclei_array)

    def test_method_calculate_returns_the_nuclear_attraction_energy_of_hydrogen_1s(self):
        kinetic_energy = self.kinetic_energy_integral.calculate(0, 0)
        testing.assert_approx_equal(kinetic_energy, -1.8804, 4)

    def test_method_calculate_returns_the_nuclear_attraction_energy_of_quantum_mechanical_off_diagonal_parts(self):
        kinetic_energy = self.kinetic_energy_integral.calculate(0, 1)
        testing.assert_approx_equal(kinetic_energy, -1.1948, 4)
