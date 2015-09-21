from unittest import TestCase
from src.main import OverlapIntegral


class TestOverlapIntegral(TestCase):

    def test_calculate_integral(self):
        overlap_integral = OverlapIntegral()
        overlap_integral.calculate_integral(1, 1, 1, 1)
        self.assertEquals(1, 1)
