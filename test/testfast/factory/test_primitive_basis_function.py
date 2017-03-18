from unittest import TestCase
from numpy import testing
from src.objects import PrimitiveBasis


class TestPrimitiveBasisFactory(TestCase):

    def test_expand_normalisation_returns_the_normalisation_of_a_gaussian_1(self):
        primitive_gaussian = PrimitiveBasis(0.15432897, 3.42525091, (0, 0, 0.7316), (0, 0, 0))
        testing.assert_approx_equal(primitive_gaussian.normalisation, 1.794441832218435, 7)

    def test_expand_normalisation_returns_the_normalisation_of_a_gaussian_2(self):
        primitive_gaussian = PrimitiveBasis(0.53532814, 0.62391373, (0, 0, 0.7316), (0, 0, 0))
        testing.assert_approx_equal(primitive_gaussian.normalisation, 0.5003264923314032, 7)

    def test_expand_normalisation_returns_the_normalisation_of_a_gaussian_3(self):
        primitive_gaussian = PrimitiveBasis(0.44463454, 0.16885540, (0, 0, 0.7316), (0, 0, 0))
        testing.assert_approx_equal(primitive_gaussian.normalisation, 0.18773545851092535, 7)

    def test_expand_normalisation_returns_the_normalisation_of_a_gaussian_4(self):
        primitive_gaussian = PrimitiveBasis(0.15432897, 9.75393461, (0, 0, 0.7316), (0, 0, 0))
        testing.assert_approx_equal(primitive_gaussian.normalisation, 3.9336432656254527, 7)

    def test_expand_normalisation_returns_the_normalisation_of_a_gaussian_5(self):
        primitive_gaussian = PrimitiveBasis(0.53532814, 1.77669115, (0, 0, 0.7316), (0, 0, 0))
        testing.assert_approx_equal(primitive_gaussian.normalisation, 1.0967787981767012, 7)

    def test_expand_normalisation_returns_the_normalisation_of_a_gaussian_6(self):
        primitive_gaussian = PrimitiveBasis(0.44463454, 0.48084429, (0, 0, 0.7316), (0, 0, 0))
        testing.assert_approx_equal(primitive_gaussian.normalisation, 0.41154131374122654, 7)
