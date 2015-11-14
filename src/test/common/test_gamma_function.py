from unittest import TestCase
from numpy import testing
from src.main.common import BoysFunction


class TestGammaFunction(TestCase):

    def test_input_of_nu_0_and_u_0_returns_1(self):
        output = BoysFunction.calculate(0, 0)
        testing.assert_approx_equal(output, 1, 6)

    def test_input_of_nu_0_and_4_802164876034457e_minus_31_returns_1(self):
        output = BoysFunction.calculate(0, 4.802164876034457e-31)
        testing.assert_approx_equal(output, 1, 6)

    def test_input_of_nu_0_and_u_4_4757221852728435_returns_0_417742(self):
        output = BoysFunction.calculate(0, 4.4757221852728435)
        testing.assert_approx_equal(output, 0.417742, 6)

    def test_input_of_nu_0_and_u_12_0769452374367_returns_0_2550(self):
        output = BoysFunction.calculate(0, 12.0769452374367)
        testing.assert_approx_equal(output, 0.255015, 6)

    def test_input_of_nu_0_and_u_130_80051761256559_returns_0_077489(self):
        output = BoysFunction.calculate(0, 130.80051761256559)
        testing.assert_approx_equal(output, 0.077489, 6)

    def test_input_of_nu_1_and_0_minus_31_returns_0_333333(self):
        output = BoysFunction.calculate(1, 0)
        testing.assert_approx_equal(output, 0.333333, 6)

    def test_input_of_nu_1_and_4_802164876034457e_minus_31_returns_0_333333(self):
        output = BoysFunction.calculate(1, 4.802164876034457e-31)
        testing.assert_approx_equal(output, 0.333333, 6)

    def test_input_of_nu_1_and_u_4_4757221852728435_returns_0_045396(self):
        output = BoysFunction.calculate(1, 4.4757221852728435)
        testing.assert_approx_equal(output, 0.045396, 6)

    def test_input_of_nu_1_and_u_12_0769452374367_returns_0_01055(self):
        output = BoysFunction.calculate(1, 12.0769452374367)
        testing.assert_approx_equal(output, 0.0105576, 6)

    def test_input_of_nu_1_and_u_130_80051761256559_returns_0_00029621(self):
        output = BoysFunction.calculate(1, 130.80051761256559)
        testing.assert_approx_equal(output, 0.000296210, 6)

    def test_input_of_nu_2_and_0_minus_31_returns_0_200000(self):
        output = BoysFunction.calculate(2, 0)
        testing.assert_approx_equal(output, 0.200000, 6)

    def test_input_of_nu_2_and_4_802164876034457e_minus_31_returns_0_200000(self):
        output = BoysFunction.calculate(2, 4.802164876034457e-31)
        testing.assert_approx_equal(output, 0.200000, 6)

    def test_input_of_nu_2_and_u_4_4757221852728435_returns_0_0139425(self):
        output = BoysFunction.calculate(2, 4.4757221852728435)
        testing.assert_approx_equal(output, 0.0139425, 6)

    def test_input_of_nu_2_and_u_12_0769452374367_returns_0_00131107(self):
        output = BoysFunction.calculate(2, 12.0769452374367)
        testing.assert_approx_equal(output, 0.00131106, 6)

    def test_input_of_nu_2_and_u_130_80051761256559_returns_3_3969e_minus_6(self):
        output = BoysFunction.calculate(2, 130.80051761256559)
        testing.assert_approx_equal(output, 3.39689e-6, 6)
