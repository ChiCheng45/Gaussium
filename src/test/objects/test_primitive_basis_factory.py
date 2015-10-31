from unittest import TestCase
from src.main.objects import PrimitiveBasisFactory
from src.main.objects import PrimitiveBasis
from numpy import testing

class TestPrimitiveBasisFactory(TestCase):


    def test_expand_basis_returns_a_array_of_coefficients_for_a_s_type_gausssian(self):
        array = PrimitiveBasisFactory().expand_basis('S', [0.15432897, 3.42525091])
        self.assertEquals(array[0].get_orbital_type(), 'S')
        self.assertEquals(array[0].get_contraction(), 0.15432897)
        self.assertEquals(array[0].get_exponent(), 3.42525091)
        self.assertEquals(array[0].get_integral_exponents(), [0, 0, 0])

    def test_del_operate_on_s_gaussian_returns_4_gaussians(self):
        primitive_gaussian = PrimitiveBasis('S', 0.15432897, 3.42525091, [0, 0, 0])
        array = PrimitiveBasisFactory().del_operator(primitive_gaussian)
        self.assertEquals(len(array), 4)

    def test_del_operate_on_gaussian_returns_s_orbital_1(self):
        primitive_gaussian = PrimitiveBasis('S', 0.15432897, 3.42525091, [0, 0, 0])
        array = PrimitiveBasisFactory().del_operator(primitive_gaussian)
        self.assertEquals(array[0].get_orbital_type(), 'S')
        testing.assert_approx_equal(array[0].get_contraction(), 1.585846, 7)
        self.assertEquals(array[0].get_exponent(), 3.42525091)
        self.assertEquals(array[0].get_integral_exponents(), [0, 0, 0])

    def test_del_operate_on_gaussian_returns_s_orbital_2(self):
        primitive_gaussian = PrimitiveBasis('S', 0.15432897, 3.42525091, [0, 0, 0])
        array = PrimitiveBasisFactory().del_operator(primitive_gaussian)
        self.assertEquals(array[1].get_orbital_type(), 'S')
        testing.assert_approx_equal(array[1].get_contraction(), -3.621281, 7)
        self.assertEquals(array[1].get_exponent(), 3.42525091)
        self.assertEquals(array[1].get_integral_exponents(), [2, 0, 0])