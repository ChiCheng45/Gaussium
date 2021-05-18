import os
from unittest import TestCase
from unittest.mock import MagicMock
from numpy import testing
from gaussium.common import read_basis_set_file


class TestFileInputBasis(TestCase):

    def setUp(self):
        helium = MagicMock(element='HELIUM', charge=2, mass=4, coordinates=(0.000000, 0.000000, 0.7316))
        hydrogen = MagicMock(element='HYDROGEN', charge=1, mass=1, coordinates=(0.000000, 0.000000, -0.7316))
        self.nuclei_array = [helium, hydrogen]

    def test_create_basis_set_array_returns_basis_array_with_coordinates_for_sto_3g(self):
        basis_file = os.path.join(os.path.dirname(__file__), 'STO-3G.gbs')
        basis_array = read_basis_set_file(basis_file, self.nuclei_array)
        testing.assert_array_equal(basis_array[0].coordinates, (0.000000, 0.000000, 0.7316))
        testing.assert_array_equal(basis_array[1].coordinates, (0.000000, 0.000000, -0.7316))

    def test_create_basis_set_array_returns_basis_array_with_s_orbital_type_exponents_for_sto_3g(self):
        basis_file = os.path.join(os.path.dirname(__file__), 'STO-3G.gbs')
        basis_array = read_basis_set_file(basis_file, self.nuclei_array)
        testing.assert_array_equal(basis_array[0].integral_exponents, (0, 0, 0))
        testing.assert_array_equal(basis_array[1].integral_exponents, (0, 0, 0))

    def test_create_basis_set_array_returns_basis_array_with_primitives_with_coordinates_for_sto_3g(self):
        basis_file = os.path.join(os.path.dirname(__file__), 'STO-3G.gbs')
        basis_array = read_basis_set_file(basis_file, self.nuclei_array)
        self.assertEqual(basis_array[0].primitive_gaussian_array[0].coordinates, (0.000000, 0.000000, 0.7316))
        self.assertEqual(basis_array[0].primitive_gaussian_array[1].coordinates, (0.000000, 0.000000, 0.7316))
        self.assertEqual(basis_array[0].primitive_gaussian_array[2].coordinates, (0.000000, 0.000000, 0.7316))
        self.assertEqual(basis_array[1].primitive_gaussian_array[0].coordinates, (0.000000, 0.000000, -0.7316))
        self.assertEqual(basis_array[1].primitive_gaussian_array[1].coordinates, (0.000000, 0.000000, -0.7316))
        self.assertEqual(basis_array[1].primitive_gaussian_array[2].coordinates, (0.000000, 0.000000, -0.7316))

    def test_create_basis_set_array_returns_basis_array_with_primitives_with_s_orbital_type_exponents_for_sto_3g(self):
        basis_file = os.path.join(os.path.dirname(__file__), 'STO-3G.gbs')
        basis_array = read_basis_set_file(basis_file, self.nuclei_array)
        self.assertEqual(basis_array[0].primitive_gaussian_array[0].integral_exponents, (0, 0, 0))
        self.assertEqual(basis_array[0].primitive_gaussian_array[1].integral_exponents, (0, 0, 0))
        self.assertEqual(basis_array[0].primitive_gaussian_array[2].integral_exponents, (0, 0, 0))
        self.assertEqual(basis_array[1].primitive_gaussian_array[0].integral_exponents, (0, 0, 0))
        self.assertEqual(basis_array[1].primitive_gaussian_array[1].integral_exponents, (0, 0, 0))
        self.assertEqual(basis_array[1].primitive_gaussian_array[2].integral_exponents, (0, 0, 0))

    def test_create_basis_set_array_returns_basis_array_with_contraction_coefficients_for_sto_3g(self):
        basis_file = os.path.join(os.path.dirname(__file__), 'STO-3G.gbs')
        basis_array = read_basis_set_file(basis_file, self.nuclei_array)
        self.assertEqual(basis_array[0].primitive_gaussian_array[0].contraction, 0.15432897)
        self.assertEqual(basis_array[0].primitive_gaussian_array[1].contraction, 0.53532814)
        self.assertEqual(basis_array[0].primitive_gaussian_array[2].contraction, 0.44463454)
        self.assertEqual(basis_array[1].primitive_gaussian_array[0].contraction, 0.15432897)
        self.assertEqual(basis_array[1].primitive_gaussian_array[1].contraction, 0.53532814)
        self.assertEqual(basis_array[1].primitive_gaussian_array[2].contraction, 0.44463454)

    def test_create_basis_set_array_returns_basis_array_with_exponents_for_sto_3g(self):
        basis_file = os.path.join(os.path.dirname(__file__), 'STO-3G.gbs')
        basis_array = read_basis_set_file(basis_file, self.nuclei_array)
        self.assertEqual(basis_array[0].primitive_gaussian_array[0].exponent, 6.36242139)
        self.assertEqual(basis_array[0].primitive_gaussian_array[1].exponent, 1.15892300)
        self.assertEqual(basis_array[0].primitive_gaussian_array[2].exponent, 0.31364979)
        self.assertEqual(basis_array[1].primitive_gaussian_array[0].exponent, 3.42525091)
        self.assertEqual(basis_array[1].primitive_gaussian_array[1].exponent, 0.62391373)
        self.assertEqual(basis_array[1].primitive_gaussian_array[2].exponent, 0.16885540)

    def test_create_basis_set_array_returns_basis_array_with_correct_coordinates_for_3_21g(self):
        basis_file = os.path.join(os.path.dirname(__file__), '3-21G.gbs')
        basis_array = read_basis_set_file(basis_file, self.nuclei_array)
        testing.assert_array_equal(basis_array[0].coordinates, (0.000000, 0.000000, 0.7316))
        testing.assert_array_equal(basis_array[1].coordinates, (0.000000, 0.000000, 0.7316))
        testing.assert_array_equal(basis_array[2].coordinates, (0.000000, 0.000000, -0.7316))
        testing.assert_array_equal(basis_array[3].coordinates, (0.000000, 0.000000, -0.7316))
