from unittest import TestCase
from unittest.mock import MagicMock
from numpy import testing
from src.main.common import Symmetry


class TestSymmetrySort(TestCase):

    def setUp(self):
        molecule = MagicMock(nuclei_array=None, point_group=None)
        self.symmetry = Symmetry(molecule, None)

    def test_sort_indices_1(self):
        out = self.symmetry.sort_index(6, 2, 1, 0)
        testing.assert_array_equal(out, (0, 1, 2, 6))

    def test_sort_indices_2(self):
        out = self.symmetry.sort_index(0, 0, 1, 0)
        testing.assert_array_equal(out, (0, 0, 0, 1))

    def test_sort_indices_3(self):
        out = self.symmetry.sort_index(1, 0, 0, 0)
        testing.assert_array_equal(out, (0, 0, 0, 1))

    def test_sort_indices_4(self):
        out = self.symmetry.sort_index(1, 0, 1, 0)
        testing.assert_array_equal(out, (0, 1, 0, 1))

    def test_sort_indices_5(self):
        out = self.symmetry.sort_index(2, 2, 1, 0)
        testing.assert_array_equal(out, (0, 1, 2, 2))

    def test_sort_indices_6(self):
        out = self.symmetry.sort_index(2, 1, 1, 0)
        testing.assert_array_equal(out, (0, 1, 1, 2))

    def test_sort_indices_7(self):
        out = self.symmetry.sort_index(1, 1, 1, 0)
        testing.assert_array_equal(out, (0, 1, 1, 1))

    def test_sort_indices_8(self):
        out = self.symmetry.sort_index(2, 1, 2, 0)
        testing.assert_array_equal(out, (0, 2, 1, 2))

    def test_sort_indices_9(self):
        out = self.symmetry.sort_index(1, 1, 0, 2)
        testing.assert_array_equal(out, (0, 2, 1, 1))
