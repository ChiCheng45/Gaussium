from unittest import TestCase
from numpy import testing
from src.main.integrals import binomial_coefficient
from src.main.integrals import combination


class TestBinomialCoefficientsFunction(TestCase):

    # returns correct answer for the coefficient of x^2 for (x + 0.7313240)(x - 0.7313240)
    def test_calculate_coefficient_return_the_correct_answer_1(self):
        coefficient = binomial_coefficient(2, 1, 1, 0.7313240, -0.7313240)
        self.assertEquals(coefficient, 1)

    # returns correct answer for the coefficient of x^1 for (x + 0.7313240)(x - 0.7313240)
    def test_calculate_coefficient_return_the_correct_answer_2(self):
        coefficient = binomial_coefficient(1, 1, 1, 0.7313240, -0.7313240)
        self.assertEquals(coefficient, 0)

    # returns correct answer for the coefficient of x^0 for (x + 0.7313240)(x - 0.7313240)
    def test_calculate_coefficient_return_the_correct_answer_3(self):
        coefficient = binomial_coefficient(0, 1, 1, 0.7313240, -0.7313240)
        testing.assert_approx_equal(coefficient, -0.534834793, 9)

    # test C(0, 0) = 1
    def test_calculate_combination_return_the_correct_answer_1(self):
        coefficient = combination(0, 0)
        self.assertEquals(coefficient, 1)

    # test C(1, 0) = 0
    def test_calculate_combination_return_the_correct_answer_2(self):
        coefficient = combination(1, 0)
        self.assertEquals(coefficient, 1)

    # test C(0, 1) = 0
    def test_calculate_combination_return_the_correct_answer_3(self):
        coefficient = combination(0, 1)
        self.assertEquals(coefficient, 0)

    # test C(1, 1) = 0
    def test_calculate_combination_return_the_correct_answer_4(self):
        coefficient = combination(1, 1)
        self.assertEquals(coefficient, 1)
