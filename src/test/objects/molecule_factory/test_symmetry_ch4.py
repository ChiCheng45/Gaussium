from unittest import TestCase
from unittest.mock import MagicMock
from numpy import testing
from src.main.objects import MoleculeFactory


class TestSymmetryCH4(TestCase):

    def setUp(self):
        carbon_1 = MagicMock(element='CARBON', charge=6, mass=12, coordinates=(-0.98781, 0.41551, 0.00000))
        hydrogen_1 = MagicMock(element='HYDROGEN', charge=1, mass=1, coordinates=(0.08219, 0.41551, 0.00000))
        hydrogen_2 = MagicMock(element='HYDROGEN', charge=1, mass=1, coordinates=(-1.34447, 0.70319, -0.96692))
        hydrogen_3 = MagicMock(element='HYDROGEN', charge=1, mass=1, coordinates=(-1.34448, 1.10904, 0.73260))
        hydrogen_4 = MagicMock(element='HYDROGEN', charge=1, mass=1, coordinates=(-1.34448, -0.56571, 0.23432))
        self.nuclei_array_ch4 = [carbon_1, hydrogen_1, hydrogen_2, hydrogen_3, hydrogen_4]
        self.molecule_factory = MoleculeFactory()

    def test_center_molecule_changes_to_carbon_1_at_origin(self):
        carbon_1 = self.molecule_factory.center_molecule(self.nuclei_array_ch4)[0]
        testing.assert_array_almost_equal(carbon_1.coordinates, (0.0, 0.0, 0.0), 6)

    def test_center_molecule_changes_to_hydrogen_1(self):
        hydrogen_1 = self.molecule_factory.center_molecule(self.nuclei_array_ch4)[1]
        testing.assert_array_almost_equal(hydrogen_1.coordinates, (1.07, 0.0, 0.0), 6)

    def test_center_molecule_changes_to_hydrogen_2(self):
        hydrogen_2 = self.molecule_factory.center_molecule(self.nuclei_array_ch4)[2]
        testing.assert_array_almost_equal(hydrogen_2.coordinates, (-0.35666, 0.28768, -0.96692), 6)

    def test_center_molecule_changes_to_hydrogen_3(self):
        hydrogen_3 = self.molecule_factory.center_molecule(self.nuclei_array_ch4)[3]
        testing.assert_array_almost_equal(hydrogen_3.coordinates, (-0.35667, 0.69353, 0.73260), 6)

    def test_center_molecule_changes_to_hydrogen_4(self):
        hydrogen_4 = self.molecule_factory.center_molecule(self.nuclei_array_ch4)[4]
        testing.assert_array_almost_equal(hydrogen_4.coordinates, (-0.35667, -0.98122, 0.23432), 6)

    def test_brute_force_rotation_symmetry_returns_list_of_seven_axis_of_rotations(self):
        nuclei_array = self.molecule_factory.center_molecule(self.nuclei_array_ch4)
        rotation, reflection = self.molecule_factory.brute_force_symmetry(nuclei_array)
        self.assertEqual(len(rotation), 7)

    def test_brute_force_rotation_symmetry_returns_list_of_four_axis_of_rotations_with_n_three(self):
        nuclei_array = self.molecule_factory.center_molecule(self.nuclei_array_ch4)
        rotation, reflection = self.molecule_factory.brute_force_symmetry(nuclei_array)
        self.assertEqual([symmetry.fold for symmetry in rotation].count(3), 4)

    def test_brute_force_rotation_symmetry_returns_list_of_four_axis_of_rotations_with_n_two(self):
        nuclei_array = self.molecule_factory.center_molecule(self.nuclei_array_ch4)
        rotation, reflection = self.molecule_factory.brute_force_symmetry(nuclei_array)
        self.assertEqual([symmetry.fold for symmetry in rotation].count(2), 3)

    def test_brute_force_reflection_symmetry_returns_a_list_of_six_reflection_planes(self):
        nuclei_array = self.molecule_factory.center_molecule(self.nuclei_array_ch4)
        rotation, reflection = self.molecule_factory.brute_force_symmetry(nuclei_array)
        self.assertEqual(len(reflection), 6)

    def test_check_linear_returns_false(self):
        nuclei_array = self.molecule_factory.center_molecule(self.nuclei_array_ch4)
        rotation, reflection = self.molecule_factory.brute_force_symmetry(nuclei_array)
        self.molecule_factory.standard_orientation(nuclei_array, rotation, reflection)
        boolean = self.molecule_factory.check_linear(nuclei_array)
        self.assertEqual(boolean, False)

    def test_check_high_symmetry_returns_true(self):
        nuclei_array = self.molecule_factory.center_molecule(self.nuclei_array_ch4)
        rotation, reflection = self.molecule_factory.brute_force_symmetry(nuclei_array)
        self.molecule_factory.standard_orientation(nuclei_array, rotation, reflection)
        boolean = self.molecule_factory.check_high_symmetry(rotation)
        self.assertEqual(boolean, True)

    def test_check_inversion_symmetry_returns_false(self):
        nuclei_array = self.molecule_factory.center_molecule(self.nuclei_array_ch4)
        rotation, reflection = self.molecule_factory.brute_force_symmetry(nuclei_array)
        self.molecule_factory.standard_orientation(nuclei_array, rotation, reflection)
        boolean = self.molecule_factory.check_inversion_symmetry(nuclei_array)
        self.assertEqual(boolean, False)

    def test_point_group_returns_t_d_symmetry_for_methane(self):
        symmetry = self.molecule_factory.point_group(self.nuclei_array_ch4).point_group
        testing.assert_equal(symmetry, 'T_{d}')