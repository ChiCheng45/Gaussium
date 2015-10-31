from unittest import TestCase
from unittest.mock import MagicMock
from src.main.common import BinomialCoefficientsFunction
from src.main.common import SFunction
import time

class TestSFunction(TestCase):

    def setUp(self):
        self.s_function = SFunction(BinomialCoefficientsFunction)
