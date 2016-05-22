from unittest import TestCase
from unittest.mock import MagicMock
from numpy import testing
from src.main.common import Symmetry


class TestSymmetrySort(TestCase):

    def test_sort_indices_1(self):
        out = Symmetry.sort(6, 2, 1, 0)
        testing.assert_array_equal(out, (0, 1, 2, 6))

    def test_sort_indices_2(self):
        out = Symmetry.sort(0, 0, 1, 0)
        testing.assert_array_equal(out, (0, 0, 0, 1))

    def test_sort_indices_3(self):
        out = Symmetry.sort(1, 0, 0, 0)
        testing.assert_array_equal(out, (0, 0, 0, 1))

    def test_sort_indices_4(self):
        out = Symmetry.sort(1, 0, 1, 0)
        testing.assert_array_equal(out, (0, 1, 0, 1))

    def test_sort_indices_5(self):
        out = Symmetry.sort(2, 2, 1, 0)
        testing.assert_array_equal(out, (0, 1, 2, 2))

    def test_sort_indices_6(self):
        out = Symmetry.sort(2, 1, 1, 0)
        testing.assert_array_equal(out, (0, 1, 1, 2))

    def test_sort_indices_7(self):
        out = Symmetry.sort(1, 1, 1, 0)
        testing.assert_array_equal(out, (0, 1, 1, 1))

    def test_sort_indices_8(self):
        out = Symmetry.sort(2, 1, 2, 0)
        testing.assert_array_equal(out, (0, 2, 1, 2))

    def test_sort_indices_9(self):
        out = Symmetry.sort(1, 1, 0, 2)
        testing.assert_array_equal(out, (0, 2, 1, 1))


class TestSymmetryChangeCoordinatesHe(TestCase):

    def setUp(self):
        helium_1 = MagicMock(element='HELIUM', charge=2, mass=4, coordinates=(-0.98781, 0.41551, 0.00000))
        self.nuclei_array_he = [helium_1]

    def test_move_nuclei_to_the_origin(self):
        helium = Symmetry.standard_orientation(self.nuclei_array_he)[0]
        testing.assert_array_equal(helium.coordinates, (0.0, 0.0, 0.0))


class TestSymmetryChangeCoordinatesN2O(TestCase):

    def setUp(self):
        nitrogen_1 = MagicMock(element='NITROGEN', charge=7, mass=14, coordinates=(0.0000000000, 0.0000000000, -2.2684205883))
        nitrogen_2 = MagicMock(element='NITROGEN', charge=7, mass=14, coordinates=(0.0000000000, 0.0000000000, -0.1349300877))
        oxygen_1 = MagicMock(element='OXYGEN', charge=8, mass=16, coordinates=(0.0000000000, 0.0000000000, 2.1042369647))
        self.nuclei_array_n2o = [nitrogen_1, nitrogen_2, oxygen_1]

    def test_change_coordinates_changes_to_standard_orientation(self):
        carbon = Symmetry.standard_orientation(self.nuclei_array_n2o)[0]
        testing.assert_array_almost_equal(carbon.coordinates, (0.0, 0.0, 0.0), 6)


class TestSymmetryChangeCoordinatesCH4(TestCase):

    def setUp(self):
        carbon_1 = MagicMock(element='CARBON', charge=6, mass=12, coordinates=(-0.98781, 0.41551, 0.00000))
        hydrogen_1 = MagicMock(element='HYDROGEN', charge=1, mass=1, coordinates=(0.08219, 0.41551, 0.00000))
        hydrogen_2 = MagicMock(element='HYDROGEN', charge=1, mass=1, coordinates=(-1.34447, 0.70319, -0.96692))
        hydrogen_3 = MagicMock(element='HYDROGEN', charge=1, mass=1, coordinates=(-1.34448, 1.10904, 0.73260))
        hydrogen_4 = MagicMock(element='HYDROGEN', charge=1, mass=1, coordinates=(-1.34448, -0.56571, 0.23432))
        self.nuclei_array_ch4 = [carbon_1, hydrogen_1, hydrogen_2, hydrogen_3, hydrogen_4]

    def test_change_coordinates_changes_to_standard_orientation(self):
       carbon = Symmetry.standard_orientation(self.nuclei_array_ch4)[0]
       testing.assert_array_almost_equal(carbon.coordinates, (0.0, 0.0, 0.0), 6)
