from unittest import TestCase
from numpy import testing
from gaussium.factory import del_operator
from gaussium.factory import expand_basis_set
from gaussium.objects import PrimitiveBasis


class TestPrimitiveBasisFactory(TestCase):

    def test_expand_basis_returns_a_array_of_coefficients_for_a_s_type_gaussian(self):
        array = expand_basis_set([['S', [0.15432897, 3.42525091]]], (0, 0, 0.7316))
        self.assertEqual(array[0].integral_exponents, (0, 0, 0))
        self.assertEqual(array[0].coordinates, (0, 0, 0.7316))
        self.assertEqual(array[0].primitive_gaussian_array[0].contraction, 0.15432897)
        self.assertEqual(array[0].primitive_gaussian_array[0].exponent, 3.42525091)

    def test_del_operate_on_s_gaussian_returns_four_gaussians(self):
        primitive_gaussian = PrimitiveBasis(0.15432897, 3.42525091, (0, 0, 0.7316), (0, 0, 0))
        array = del_operator(primitive_gaussian)
        self.assertEqual(len(array), 4)

    def test_del_operate_on_gaussian_returns_s_orbital_1(self):
        primitive_gaussian = PrimitiveBasis(0.15432897, 3.42525091, (0, 0, 0.7316), (0, 0, 0))
        array = del_operator(primitive_gaussian)
        testing.assert_approx_equal(array[0].contraction, 10.27575273, 7)
        self.assertEqual(array[0].exponent, 3.42525091)
        self.assertEqual(array[0].integral_exponents, (0, 0, 0))
        self.assertEqual(array[1].coordinates, (0, 0, 0.7316))

    def test_del_operate_on_gaussian_returns_s_orbital_2(self):
        primitive_gaussian = PrimitiveBasis(0.15432897, 3.42525091, (0, 0, 0.7316), (0, 0, 0))
        array = del_operator(primitive_gaussian)
        testing.assert_approx_equal(array[1].contraction, -23.46468759, 7)
        self.assertEqual(array[1].exponent, 3.42525091)
        self.assertEqual(array[1].integral_exponents, (2, 0, 0))
        self.assertEqual(array[1].coordinates, (0, 0, 0.7316))
