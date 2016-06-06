from unittest import TestCase
from unittest.mock import MagicMock
from numpy import testing
from src.main.objects import MoleculeFactory


class TestSymmetryChangeCoordinatesHe(TestCase):

    def setUp(self):
        helium_1 = MagicMock(element='HELIUM', charge=2, mass=4, coordinates=(-0.98781, 0.41551, 0.00000))
        self.nuclei_array_he = [helium_1]

    def test_move_nuclei_to_the_origin(self):
        helium = MoleculeFactory.point_group(self.nuclei_array_he)[0]
        testing.assert_array_equal(helium.coordinates, (0.0, 0.0, 0.0))


class TestSymmetryChangeCoordinatesN2O(TestCase):

    def setUp(self):
        nitrogen_1 = MagicMock(element='NITROGEN', charge=7, mass=14, coordinates=(0.0000000000, 0.0000000000, -2.2684205883))
        nitrogen_2 = MagicMock(element='NITROGEN', charge=7, mass=14, coordinates=(0.0000000000, 0.0000000000, -0.1349300877))
        oxygen_1 = MagicMock(element='OXYGEN', charge=8, mass=16, coordinates=(0.0000000000, 0.0000000000, 2.1042369647))
        self.nuclei_array_n2o = [nitrogen_1, nitrogen_2, oxygen_1]

    # def test_change_coordinates_changes_to_standard_orientation(self):
    #     nitrogen_1 = MoleculeFactory.point_group(self.nuclei_array_n2o).nuclei_array[0]
    #     testing.assert_array_almost_equal(nitrogen_1.coordinates, (0.0, 0.0, 0.0), 6)


class TestSymmetryChangeCoordinatesCH4(TestCase):

    def setUp(self):
        carbon_1 = MagicMock(element='CARBON', charge=6, mass=12, coordinates=(-0.98781, 0.41551, 0.00000))
        hydrogen_1 = MagicMock(element='HYDROGEN', charge=1, mass=1, coordinates=(0.08219, 0.41551, 0.00000))
        hydrogen_2 = MagicMock(element='HYDROGEN', charge=1, mass=1, coordinates=(-1.34447, 0.70319, -0.96692))
        hydrogen_3 = MagicMock(element='HYDROGEN', charge=1, mass=1, coordinates=(-1.34448, 1.10904, 0.73260))
        hydrogen_4 = MagicMock(element='HYDROGEN', charge=1, mass=1, coordinates=(-1.34448, -0.56571, 0.23432))
        self.nuclei_array_ch4 = [carbon_1, hydrogen_1, hydrogen_2, hydrogen_3, hydrogen_4]

    def test_center_molecule_changes_to_carbon_1_at_origin(self):
        carbon_1 = MoleculeFactory.center_molecule(self.nuclei_array_ch4)[0]
        testing.assert_array_almost_equal(carbon_1.coordinates, (0.0, 0.0, 0.0), 6)

    def test_center_molecule_changes_to_hydrogen_1(self):
        hydrogen_1 = MoleculeFactory.center_molecule(self.nuclei_array_ch4)[1]
        testing.assert_array_almost_equal(hydrogen_1.coordinates, (1.07, 0.0, 0.0), 6)

    def test_center_molecule_changes_to_hydrogen_2(self):
        hydrogen_2 = MoleculeFactory.center_molecule(self.nuclei_array_ch4)[2]
        testing.assert_array_almost_equal(hydrogen_2.coordinates, (-0.35666, 0.28768, -0.96692), 6)

    def test_center_molecule_changes_to_hydrogen_3(self):
        hydrogen_3 = MoleculeFactory.center_molecule(self.nuclei_array_ch4)[3]
        testing.assert_array_almost_equal(hydrogen_3.coordinates, (-0.35667, 0.69353, 0.73260), 6)

    def test_center_molecule_changes_to_hydrogen_4(self):
        hydrogen_4 = MoleculeFactory.center_molecule(self.nuclei_array_ch4)[4]
        testing.assert_array_almost_equal(hydrogen_4.coordinates, (-0.35667, -0.98122, 0.23432), 6)

    def test_brute_force_rotation_symmetry_returns_list_of_seven_axis_of_rotations(self):
        nuclei_array = MoleculeFactory.center_molecule(self.nuclei_array_ch4)
        symmetry_list = MoleculeFactory.brute_force_rotation_symmetry(nuclei_array)
        self.assertEqual(len(symmetry_list), 7)

    def test_point_group_returns_t_d_symmetry_for_methane(self):
        symmetry = MoleculeFactory.point_group(self.nuclei_array_ch4).symmetry
        testing.assert_equal(symmetry, 'T_{d}')
