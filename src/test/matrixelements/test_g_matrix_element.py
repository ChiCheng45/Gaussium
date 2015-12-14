from unittest import TestCase
from src.main.matrixelements import GMatrixElement
from numpy import testing


class TestGMatrixElements(TestCase):

    def setUp(self):
        pass

    def test_sort_indices_1(self):
        out = GMatrixElement.sort(6, 2, 1, 0)
        testing.assert_array_equal(out, (0, 1, 2, 6))

    def test_sort_indices_2(self):
        out = GMatrixElement.sort(0, 0, 1, 0)
        testing.assert_array_equal(out, (0, 0, 0, 1))

    def test_sort_indices_3(self):
        out = GMatrixElement.sort(1, 0, 0, 0)
        testing.assert_array_equal(out, (0, 0, 0, 1))

    def test_sort_indices_4(self):
        out = GMatrixElement.sort(1, 0, 1, 0)
        testing.assert_array_equal(out, (0, 1, 0, 1))

    def test_sort_indices_5(self):
        out = GMatrixElement.sort(2, 2, 1, 0)
        testing.assert_array_equal(out, (0, 1, 2, 2))

    def test_sort_indices_6(self):
        out = GMatrixElement.sort(2, 1, 1, 0)
        testing.assert_array_equal(out, (0, 1, 1, 2))

    def test_sort_indices_7(self):
        out = GMatrixElement.sort(1, 1, 1, 0)
        testing.assert_array_equal(out, (0, 1, 1, 1))

    def test_sort_indices_8(self):
        out = GMatrixElement.sort(2, 1, 2, 0)
        testing.assert_array_equal(out, (0, 2, 1, 2))

    def test_sort_indices_9(self):
        out = GMatrixElement.sort(1, 1, 0, 2)
        testing.assert_array_equal(out, (0, 2, 1, 1))
