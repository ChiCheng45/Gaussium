from unittest import TestCase
from unittest.mock import MagicMock
from numpy import testing
from src.main.objects.molecule_factory import *


class TestSymmetryHe(TestCase):
    def setUp(self):
        helium_1 = MagicMock(element='HELIUM', charge=2, mass=4, coordinates=(-0.98781, 0.41551, 0.00000))
        self.nuclei_array_he = [helium_1]

    def test_move_nuclei_to_the_origin(self):
        helium = point_group(self.nuclei_array_he).nuclei_array[0]
        testing.assert_array_equal(helium.coordinates, (0.0, 0.0, 0.0))

    def test_point_group_returns_c_1_symmetry_for_helium(self):
        symmetry = point_group(self.nuclei_array_he).symmetry_group
        testing.assert_equal(symmetry, 'C_{1}')


class TestSymmetryHOF(TestCase):
    def setUp(self):
        oxygen_1 = MagicMock(element='OXYGEN', charge=8, mass=16, coordinates=(-1.4186923158, 0.1090030362, 0.0000000000))
        hydrogen_1 = MagicMock(element='HYDROGEN', charge=1, mass=1, coordinates=(-1.7313653816, -1.6895740638, 0.0000000000))
        fluorine_1 = MagicMock(element='FLUORINE', charge=9, mass=19, coordinates=(1.2899273141, 0.0031592817, 0.0000000000))
        self.nuclei_array_hof = [oxygen_1, hydrogen_1, fluorine_1]

    def test_brute_force_rotation_symmetry_returns_list_of_zero_axis_of_rotations(self):
        nuclei_array = center_molecule(self.nuclei_array_hof)
        rotation, reflection = brute_force_symmetry(nuclei_array)
        self.assertEqual(len(rotation), 0)

    def test_brute_force_reflection_symmetry_returns_list_of_one_reflection_planes(self):
        nuclei_array = center_molecule(self.nuclei_array_hof)
        rotation, reflection = brute_force_symmetry(nuclei_array)
        self.assertEqual(len(reflection), 1)

    def test_check_linear_returns_false(self):
        nuclei_array = center_molecule(self.nuclei_array_hof)
        rotation, reflection = brute_force_symmetry(nuclei_array)
        nuclei_array, rotation, reflection = standard_orientation(nuclei_array, rotation, reflection)
        boolean = check_linear(nuclei_array)
        self.assertEqual(boolean, False)

    def test_check_high_symmetry_returns_false(self):
        nuclei_array = center_molecule(self.nuclei_array_hof)
        rotation, reflection = brute_force_symmetry(nuclei_array)
        nuclei_array, rotation, reflection = standard_orientation(nuclei_array, rotation, reflection)
        boolean = check_high_symmetry(rotation)
        self.assertEqual(boolean, False)

    def test_check_sigma_h_returns_true(self):
        nuclei_array = center_molecule(self.nuclei_array_hof)
        rotation, reflection = brute_force_symmetry(nuclei_array)
        nuclei_array, rotation, reflection = standard_orientation(nuclei_array, rotation, reflection)
        boolean = check_sigma_h(reflection)
        self.assertEqual(boolean, True)

    def test_point_group_returns_c_s_symmetry_for_hypofluorous_acid(self):
        symmetry = point_group(self.nuclei_array_hof).symmetry_group
        testing.assert_equal(symmetry, 'C_{s}')


class TestSymmetryH2O(TestCase):
    def setUp(self):
        oxygen_1 = MagicMock(element='OXYGEN', charge=8, mass=16, coordinates=(0.0000000000, 0.0000000000, -0.1363928482))
        hydrogen_1 = MagicMock(element='HYDROGEN', charge=1, mass=1, coordinates=(0.0000000000, 1.4236595095, 0.9813433754))
        hydrogen_2 = MagicMock(element='HYDROGEN', charge=1, mass=1, coordinates=(0.0000000000, -1.4236595095, 0.9813433754))
        self.nuclei_array_h2o = [oxygen_1, hydrogen_1, hydrogen_2]

    def test_brute_force_rotation_symmetry_returns_list_of_one_axis_of_rotations(self):
        nuclei_array = center_molecule(self.nuclei_array_h2o)
        rotation, reflection = brute_force_symmetry(nuclei_array)
        self.assertEqual(len(rotation), 1)

    def test_brute_force_rotation_symmetry_returns_axis_of_rotation_of_two_fold_symmetry(self):
        nuclei_array = center_molecule(self.nuclei_array_h2o)
        rotation, reflection = brute_force_symmetry(nuclei_array)
        self.assertEqual(rotation[0].fold, 2)

    def test_brute_force_reflection_symmetry_returns_list_of_two_reflection_planes(self):
        nuclei_array = center_molecule(self.nuclei_array_h2o)
        rotation, reflection = brute_force_symmetry(nuclei_array)
        self.assertEqual(len(reflection), 2)

    def test_check_linear_returns_false(self):
        nuclei_array = center_molecule(self.nuclei_array_h2o)
        rotation, reflection = brute_force_symmetry(nuclei_array)
        nuclei_array, rotation, reflection = standard_orientation(nuclei_array, rotation, reflection)
        boolean = check_linear(nuclei_array)
        self.assertEqual(boolean, False)

    def test_check_high_symmetry_returns_false(self):
        nuclei_array = center_molecule(self.nuclei_array_h2o)
        rotation, reflection = brute_force_symmetry(nuclei_array)
        nuclei_array, rotation, reflection = standard_orientation(nuclei_array, rotation, reflection)
        boolean = check_high_symmetry(rotation)
        self.assertEqual(boolean, False)

    def test_get_n_fold_returns_two(self):
        nuclei_array = center_molecule(self.nuclei_array_h2o)
        rotation, reflection = brute_force_symmetry(nuclei_array)
        nuclei_array, rotation, reflection = standard_orientation(nuclei_array, rotation, reflection)
        n = get_n_fold(rotation)
        self.assertEqual(n, 2)

    def test_check_n_two_fold_rotation_perpendicular_to_n_fold_returns_false(self):
        nuclei_array = center_molecule(self.nuclei_array_h2o)
        rotation, reflection = brute_force_symmetry(nuclei_array)
        nuclei_array, rotation, reflection = standard_orientation(nuclei_array, rotation, reflection)
        boolean = check_n_two_fold_perpendicular_to_n_fold(rotation)
        self.assertEqual(boolean, False)

    def test_check_sigma_h_returns_false(self):
        nuclei_array = center_molecule(self.nuclei_array_h2o)
        rotation, reflection = brute_force_symmetry(nuclei_array)
        nuclei_array, rotation, reflection = standard_orientation(nuclei_array, rotation, reflection)
        boolean = check_sigma_h(reflection)
        self.assertEqual(boolean, False)

    def test_check_n_sigma_v_returns_true(self):
        nuclei_array = center_molecule(self.nuclei_array_h2o)
        rotation, reflection = brute_force_symmetry(nuclei_array)
        nuclei_array, rotation, reflection = standard_orientation(nuclei_array, rotation, reflection)
        boolean = check_n_sigma_v(2, reflection)
        self.assertEqual(boolean, True)

    def test_point_group_returns_c_2v_symmetry_for_water(self):
        symmetry = point_group(self.nuclei_array_h2o).symmetry_group
        testing.assert_equal(symmetry, 'C_{2v}')


class TestSymmetryNH3(TestCase):
    pass


class TestSymmetryC2H4(TestCase):
    def setUp(self):
        carbon_1 = MagicMock(element='CARBON', charge=6, mass=12, coordinates=(0.0000000000, 1.2594652672, 0.0000000000))
        carbon_2 = MagicMock(element='CARBON', charge=6, mass=12, coordinates=(0.0000000000, -1.2594652672, 0.0000000000))
        hydrogen_1 = MagicMock(element='HYDROGEN', charge=1, mass=1, coordinates=(1.7400646600, 2.3216269636, 0.0000000000))
        hydrogen_2 = MagicMock(element='HYDROGEN', charge=1, mass=1, coordinates=(-1.7400646600, 2.3216269636, 0.0000000000))
        hydrogen_3 = MagicMock(element='HYDROGEN', charge=1, mass=1, coordinates=(1.7400646600, -2.3216269636, 0.0000000000))
        hydrogen_4 = MagicMock(element='HYDROGEN', charge=1, mass=1, coordinates=(-1.7400646600, -2.3216269636, 0.0000000000))
        self.nuclei_array_c2h4 = [carbon_1, carbon_2, hydrogen_1, hydrogen_2, hydrogen_3, hydrogen_4]

    def test_brute_force_rotation_symmetry_returns_list_of_three_axis_of_rotations(self):
        nuclei_array = center_molecule(self.nuclei_array_c2h4)
        rotation, reflection = brute_force_symmetry(nuclei_array)
        self.assertEqual(len(rotation), 3)

    def test_brute_force_reflection_symmetry_returns_list_of_three_reflection_planes(self):
        nuclei_array = center_molecule(self.nuclei_array_c2h4)
        rotation, reflection = brute_force_symmetry(nuclei_array)
        self.assertEqual(len(reflection), 3)

    def test_check_inversion_symmetry_returns_true(self):
        nuclei_array = center_molecule(self.nuclei_array_c2h4)
        rotation, reflection = brute_force_symmetry(nuclei_array)
        nuclei_array, rotation, reflection = standard_orientation(nuclei_array, rotation, reflection)
        boolean = check_inversion_symmetry(nuclei_array)
        self.assertEqual(boolean, True)

    def test_check_linear_returns_false(self):
        nuclei_array = center_molecule(self.nuclei_array_c2h4)
        rotation, reflection = brute_force_symmetry(nuclei_array)
        nuclei_array, rotation, reflection = standard_orientation(nuclei_array, rotation, reflection)
        boolean = check_linear(nuclei_array)
        self.assertEqual(boolean, False)

    def test_check_high_symmetry_returns_false(self):
        nuclei_array = center_molecule(self.nuclei_array_c2h4)
        rotation, reflection = brute_force_symmetry(nuclei_array)
        nuclei_array, rotation, reflection = standard_orientation(nuclei_array, rotation, reflection)
        boolean = check_high_symmetry(rotation)
        self.assertEqual(boolean, False)

    def test_check_n_two_fold_rotation_perpendicular_to_n_fold_returns_true(self):
        nuclei_array = center_molecule(self.nuclei_array_c2h4)
        rotation, reflection = brute_force_symmetry(nuclei_array)
        nuclei_array, rotation, reflection = standard_orientation(nuclei_array, rotation, reflection)
        boolean = check_n_two_fold_perpendicular_to_n_fold(rotation)
        self.assertEqual(boolean, True)

    def test_check_sigma_h_returns_true(self):
        nuclei_array = center_molecule(self.nuclei_array_c2h4)
        rotation, reflection = brute_force_symmetry(nuclei_array)
        nuclei_array, rotation, reflection = standard_orientation(nuclei_array, rotation, reflection)
        boolean = check_sigma_h(reflection)
        self.assertEqual(boolean, True)

    def test_point_group_returns_d_2h_symmetry_for_ethene(self):
        symmetry = point_group(self.nuclei_array_c2h4).symmetry_group
        testing.assert_equal(symmetry, 'D_{2h}')


class TestSymmetryN2O(TestCase):
    def setUp(self):
        nitrogen_1 = MagicMock(element='NITROGEN', charge=7, mass=14, coordinates=(0.0000000000, 0.0000000000, -2.2684205883))
        nitrogen_2 = MagicMock(element='NITROGEN', charge=7, mass=14, coordinates=(0.0000000000, 0.0000000000, -0.1349300877))
        oxygen_1 = MagicMock(element='OXYGEN', charge=8, mass=16, coordinates=(0.0000000000, 0.0000000000, 2.1042369647))
        self.nuclei_array_n2o = [nitrogen_1, nitrogen_2, oxygen_1]

    def test_check_linear_returns_true(self):
        nuclei_array = center_molecule(self.nuclei_array_n2o)
        rotation, reflection = brute_force_symmetry(nuclei_array)
        nuclei_array, rotation, reflection = standard_orientation(nuclei_array, rotation, reflection)
        boolean = check_linear(nuclei_array)
        self.assertEqual(boolean, True)

    def test_check_inversion_symmetry_returns_true(self):
        nuclei_array = center_molecule(self.nuclei_array_n2o)
        rotation, reflection = brute_force_symmetry(nuclei_array)
        nuclei_array, rotation, reflection = standard_orientation(nuclei_array, rotation, reflection)
        boolean = check_inversion_symmetry(nuclei_array)
        self.assertEqual(boolean, False)

    def test_point_group_returns_c_inf_v_symmetry_for_nitrous_oxide(self):
        symmetry = point_group(self.nuclei_array_n2o).symmetry_group
        testing.assert_equal(symmetry, 'C_{inf v}')


class TestSymmetryN2(TestCase):
    def setUp(self):
        nitrogen_1 = MagicMock(element='NITROGEN', charge=7, mass=14, coordinates=(0.0000000000, 0.0000000000, 1.0399092291))
        nitrogen_2 = MagicMock(element='NITROGEN', charge=7, mass=14, coordinates=(0.0000000000, 0.0000000000, -1.0399092291))
        self.nuclei_array_n2 = [nitrogen_1, nitrogen_2]

    def test_check_linear_returns_true(self):
        nuclei_array = center_molecule(self.nuclei_array_n2)
        rotation, reflection = brute_force_symmetry(nuclei_array)
        nuclei_array, rotation, reflection = standard_orientation(nuclei_array, rotation, reflection)
        boolean = check_linear(nuclei_array)
        self.assertEqual(boolean, True)

    def test_check_inversion_symmetry_returns_true(self):
        nuclei_array = center_molecule(self.nuclei_array_n2)
        rotation, reflection = brute_force_symmetry(nuclei_array)
        nuclei_array, rotation, reflection = standard_orientation(nuclei_array, rotation, reflection)
        boolean = check_inversion_symmetry(nuclei_array)
        self.assertEqual(boolean, True)

    def test_point_group_returns_d_inf_h_symmetry_for_nitrogen(self):
        symmetry = point_group(self.nuclei_array_n2).symmetry_group
        testing.assert_equal(symmetry, 'D_{inf h}')


class TestSymmetryCH4(TestCase):
    def setUp(self):
        carbon_1 = MagicMock(element='CARBON', charge=6, mass=12, coordinates=(-0.98781, 0.41551, 0.00000))
        hydrogen_1 = MagicMock(element='HYDROGEN', charge=1, mass=1, coordinates=(0.08219, 0.41551, 0.00000))
        hydrogen_2 = MagicMock(element='HYDROGEN', charge=1, mass=1, coordinates=(-1.34447, 0.70319, -0.96692))
        hydrogen_3 = MagicMock(element='HYDROGEN', charge=1, mass=1, coordinates=(-1.34448, 1.10904, 0.73260))
        hydrogen_4 = MagicMock(element='HYDROGEN', charge=1, mass=1, coordinates=(-1.34448, -0.56571, 0.23432))
        self.nuclei_array_ch4 = [carbon_1, hydrogen_1, hydrogen_2, hydrogen_3, hydrogen_4]

    def test_center_molecule_changes_to_carbon_1_at_origin(self):
        carbon_1 = center_molecule(self.nuclei_array_ch4)[0]
        testing.assert_array_almost_equal(carbon_1.coordinates, (0.0, 0.0, 0.0), 6)

    def test_center_molecule_changes_to_hydrogen_1(self):
        hydrogen_1 = center_molecule(self.nuclei_array_ch4)[1]
        testing.assert_array_almost_equal(hydrogen_1.coordinates, (1.07, 0.0, 0.0), 6)

    def test_center_molecule_changes_to_hydrogen_2(self):
        hydrogen_2 = center_molecule(self.nuclei_array_ch4)[2]
        testing.assert_array_almost_equal(hydrogen_2.coordinates, (-0.35666, 0.28768, -0.96692), 6)

    def test_center_molecule_changes_to_hydrogen_3(self):
        hydrogen_3 = center_molecule(self.nuclei_array_ch4)[3]
        testing.assert_array_almost_equal(hydrogen_3.coordinates, (-0.35667, 0.69353, 0.73260), 6)

    def test_center_molecule_changes_to_hydrogen_4(self):
        hydrogen_4 = center_molecule(self.nuclei_array_ch4)[4]
        testing.assert_array_almost_equal(hydrogen_4.coordinates, (-0.35667, -0.98122, 0.23432), 6)

    def test_brute_force_rotation_symmetry_returns_list_of_seven_axis_of_rotations(self):
        nuclei_array = center_molecule(self.nuclei_array_ch4)
        rotation, reflection = brute_force_symmetry(nuclei_array)
        self.assertEqual(len(rotation), 7)

    def test_brute_force_rotation_symmetry_returns_list_of_four_axis_of_rotations_with_n_three(self):
        nuclei_array = center_molecule(self.nuclei_array_ch4)
        rotation, reflection = brute_force_symmetry(nuclei_array)
        self.assertEqual([symmetry.fold for symmetry in rotation].count(3), 4)

    def test_brute_force_rotation_symmetry_returns_list_of_four_axis_of_rotations_with_n_two(self):
        nuclei_array = center_molecule(self.nuclei_array_ch4)
        rotation, reflection = brute_force_symmetry(nuclei_array)
        self.assertEqual([symmetry.fold for symmetry in rotation].count(2), 3)

    def test_brute_force_reflection_symmetry_returns_a_list_of_six_reflection_planes(self):
        nuclei_array = center_molecule(self.nuclei_array_ch4)
        rotation, reflection = brute_force_symmetry(nuclei_array)
        self.assertEqual(len(reflection), 6)

    def test_check_linear_returns_false(self):
        nuclei_array = center_molecule(self.nuclei_array_ch4)
        rotation, reflection = brute_force_symmetry(nuclei_array)
        nuclei_array, rotation, reflection = standard_orientation(nuclei_array, rotation, reflection)
        boolean = check_linear(nuclei_array)
        self.assertEqual(boolean, False)

    def test_check_high_symmetry_returns_true(self):
        nuclei_array = center_molecule(self.nuclei_array_ch4)
        rotation, reflection = brute_force_symmetry(nuclei_array)
        nuclei_array, rotation, reflection = standard_orientation(nuclei_array, rotation, reflection)
        boolean = check_high_symmetry(rotation)
        self.assertEqual(boolean, True)

    def test_check_inversion_symmetry_returns_false(self):
        nuclei_array = center_molecule(self.nuclei_array_ch4)
        rotation, reflection = brute_force_symmetry(nuclei_array)
        nuclei_array, rotation, reflection = standard_orientation(nuclei_array, rotation, reflection)
        boolean = check_inversion_symmetry(nuclei_array)
        self.assertEqual(boolean, False)

    def test_point_group_returns_t_d_symmetry_for_methane(self):
        symmetry = point_group(self.nuclei_array_ch4).symmetry_group
        testing.assert_equal(symmetry, 'T_{d}')


class TestSymmetryC8H8(TestCase):
    def setUp(self):
        carbon_1 = MagicMock(element='CARBON', charge=6, mass=12, coordinates=(1.2455, 0.5367, -0.0729))
        carbon_2 = MagicMock(element='CARBON', charge=6, mass=12, coordinates=(0.9239, -0.9952, 0.0237))
        carbon_3 = MagicMock(element='CARBON', charge=6, mass=12, coordinates=(-0.1226, -0.7041, 1.1548))
        carbon_4 = MagicMock(element='CARBON', charge=6, mass=12, coordinates=(0.1989, 0.8277, 1.0582))
        carbon_5 = MagicMock(element='CARBON', charge=6, mass=12, coordinates=(0.1226, 0.7042, -1.1548))
        carbon_6 = MagicMock(element='CARBON', charge=6, mass=12, coordinates=(-0.9239, 0.9952, -0.0237))
        carbon_7 = MagicMock(element='CARBON', charge=6, mass=12, coordinates=(-1.2454, -0.5367, 0.0729))
        carbon_8 = MagicMock(element='CARBON', charge=6, mass=12, coordinates=(-0.1989, -0.8277, -1.0582))
        hydrogen_1 = MagicMock(element='HYDROGEN', charge=1, mass=1, coordinates=(2.2431, 0.9666, -0.1313))
        hydrogen_2 = MagicMock(element='HYDROGEN', charge=1, mass=1, coordinates=(1.6638, -1.7924, 0.0426))
        hydrogen_3 = MagicMock(element='HYDROGEN', charge=1, mass=1, coordinates=(-0.2209, -1.2683, 2.0797))
        hydrogen_4 = MagicMock(element='HYDROGEN', charge=1, mass=1, coordinates=(0.3583, 1.4907, 1.9059))
        hydrogen_5 = MagicMock(element='HYDROGEN', charge=1, mass=1, coordinates=(0.2208, 1.2681, -2.0799))
        hydrogen_6 = MagicMock(element='HYDROGEN', charge=1, mass=1, coordinates=(-1.6640, 1.7922, -0.0427))
        hydrogen_7 = MagicMock(element='HYDROGEN', charge=1, mass=1, coordinates=(-2.2430, -0.9665, 0.1313))
        hydrogen_8 = MagicMock(element='HYDROGEN', charge=1, mass=1, coordinates=(-0.3583, -1.4906, -1.9058))

        self.nuclei_array_c8h8 = [carbon_1, carbon_2, carbon_3, carbon_4, carbon_5, carbon_6, carbon_7, carbon_8,
        hydrogen_1, hydrogen_2, hydrogen_3, hydrogen_4, hydrogen_5, hydrogen_6, hydrogen_7, hydrogen_8]

    def test_check_inversion_symmetry_returns_true(self):
        nuclei_array = center_molecule(self.nuclei_array_c8h8)
        rotation, reflection = brute_force_symmetry(nuclei_array)
        nuclei_array, rotation, reflection = standard_orientation(nuclei_array, rotation, reflection)
        boolean = check_inversion_symmetry(nuclei_array)
        self.assertEqual(boolean, True)

    def test_point_group_returns_i_h_symmetry_for_cubane(self):
        symmetry = point_group(self.nuclei_array_c8h8).symmetry_group
        testing.assert_equal(symmetry, 'O_{h}')
