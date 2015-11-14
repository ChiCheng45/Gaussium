from unittest import TestCase
from src.main.matrixelements import GMatrixElement
from numpy import testing
import numpy as np


class TestGMatrixElements(TestCase):

    def setUp(self):
        pass

    def test_sort1(self):
        inp = (6, 2, 1, 0)
        out = GMatrixElement.sort(inp)
        testing.assert_array_equal(out, (1, 0, 6, 2))

    def test_sort2(self):
        inp = (0, 0, 1, 0)
        out = GMatrixElement.sort(inp)
        testing.assert_array_equal(out, (0, 0, 0, 1))

    def test_sort4(self):
        inp = (1, 0, 0, 0)
        out = GMatrixElement.sort(inp)
        testing.assert_array_equal(out, (0, 0, 0, 1))

    def test_sort3(self):
        inp = (1, 0, 1, 0)
        out = GMatrixElement.sort(inp)
        testing.assert_array_equal(out, (0, 1, 0, 1))

    def test_sort5(self):
        inp = (2, 2, 1, 0)
        out = GMatrixElement.sort(inp)
        testing.assert_array_equal(out, (0, 1, 2, 2))

    def test_sort6(self):
        inp = (2, 1, 1, 0)
        out = GMatrixElement.sort(inp)
        testing.assert_array_equal(out, (0, 1, 1, 2))

    def test_sort7(self):
        inp = (1, 1, 1, 0)
        out = GMatrixElement.sort(inp)
        testing.assert_array_equal(out, (0, 1, 1, 1))

    def test_sort8(self):
        inp = (2, 1, 2, 0)
        out = GMatrixElement.sort(inp)
        testing.assert_array_equal(out, (0, 2, 1, 2))

    def test_sort9(self):
        inp = (1, 1, 0, 2)
        out = GMatrixElement.sort(inp)
        testing.assert_array_equal(out, (0, 2, 1, 1))