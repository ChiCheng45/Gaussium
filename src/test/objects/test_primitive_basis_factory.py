from unittest import TestCase
from src.main.objects import PrimitiveBasisFactory

class TestPrimitiveBasisFactory(TestCase):

    def setUp(self):
        pass

    def test_expand_basis_returns_a_array_of_coefficients_for_a_s_type_gausssian(self):
        array = PrimitiveBasisFactory().expand_basis('S', [0.15432897, 3.42525091])
        self.assertEquals(array[0].get_orbital_type(), 'S')
        self.assertEquals(array[0].get_contraction(), 0.15432897)
        self.assertEquals(array[0].get_exponent(), 3.42525091)
        self.assertEquals(array[0].get_integral_exponents(), [0, 0, 0])
