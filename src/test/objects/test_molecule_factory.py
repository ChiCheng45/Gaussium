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


class TestSymmetryHOF(TestCase):

    def setUp(self):
        oxygen_1 = MagicMock(element='OXYGEN', charge=8, mass=16, coordinates=(-1.4186923158, 0.1090030362, 0.0000000000))
        hydrogen_1 = MagicMock(element='HYDROGEN', charge=1, mass=1, coordinates=(-1.7313653816, -1.6895740638, 0.0000000000))
        fluorine_1 = MagicMock(element='FLUORINE', charge=9, mass=19, coordinates=(1.2899273141, 0.0031592817, 0.0000000000))
        self.nuclei_array_hof = [oxygen_1, hydrogen_1, fluorine_1]

    def test_brute_force_rotation_symmetry_returns_list_of_zero_axis_of_rotations(self):
        nuclei_array = MoleculeFactory.center_molecule(self.nuclei_array_hof)
        symmetry_list = MoleculeFactory.brute_force_rotation_symmetry(nuclei_array)
        self.assertEqual(len(symmetry_list), 0)

    def test_brute_force_reflection_symmetry_returns_list_of_one_reflection_planes(self):
        nuclei_array = MoleculeFactory.center_molecule(self.nuclei_array_hof)
        rotation_symmetry_list = MoleculeFactory.brute_force_rotation_symmetry(nuclei_array)
        reflection_symmetry_list = MoleculeFactory.brute_force_reflection_symmetry(nuclei_array, rotation_symmetry_list)
        self.assertEqual(len(reflection_symmetry_list), 1)

    def test_check_linear_returns_false(self):
        nuclei_array = MoleculeFactory.center_molecule(self.nuclei_array_hof)
        boolean = MoleculeFactory.check_linear(nuclei_array)
        self.assertEqual(boolean, False)

    def test_check_high_symmetry_returns_false(self):
        nuclei_array = MoleculeFactory.center_molecule(self.nuclei_array_hof)
        symmetry_list = MoleculeFactory.brute_force_rotation_symmetry(nuclei_array)
        boolean = MoleculeFactory.check_high_symmetry(symmetry_list)
        self.assertEqual(boolean, False)

    def test_check_sigma_h_returns_true(self):
        nuclei_array = MoleculeFactory.center_molecule(self.nuclei_array_hof)
        rotation_symmetry_list = MoleculeFactory.brute_force_rotation_symmetry(nuclei_array)
        reflection_symmetry_list = MoleculeFactory.brute_force_reflection_symmetry(nuclei_array, rotation_symmetry_list)
        boolean = MoleculeFactory.check_sigma_h(reflection_symmetry_list)
        self.assertEqual(boolean, True)

    def test_point_group_returns_c_s_symmetry_for_hypofluorous_acid(self):
        symmetry = MoleculeFactory.point_group(self.nuclei_array_hof).symmetry_group
        testing.assert_equal(symmetry, 'C_{s}')


class TestSymmetryH2O(TestCase):

    def setUp(self):
        oxygen_1 = MagicMock(element='OXYGEN', charge=8, mass=16, coordinates=(0.0000000000, 0.0000000000, -0.1363928482))
        hydrogen_1 = MagicMock(element='HYDROGEN', charge=1, mass=1, coordinates=(0.0000000000, 1.4236595095, 0.9813433754))
        hydrogen_2 = MagicMock(element='HYDROGEN', charge=1, mass=1, coordinates=(0.0000000000, -1.4236595095, 0.9813433754))
        self.nuclei_array_h2o = [oxygen_1, hydrogen_1, hydrogen_2]

    def test_brute_force_rotation_symmetry_returns_list_of_one_axis_of_rotations(self):
        nuclei_array = MoleculeFactory.center_molecule(self.nuclei_array_h2o)
        symmetry_list = MoleculeFactory.brute_force_rotation_symmetry(nuclei_array)
        self.assertEqual(len(symmetry_list), 1)

    def test_brute_force_rotation_symmetry_returns_axis_of_rotation_of_two_fold_symmetry(self):
        nuclei_array = MoleculeFactory.center_molecule(self.nuclei_array_h2o)
        symmetry_list = MoleculeFactory.brute_force_rotation_symmetry(nuclei_array)
        self.assertEqual(symmetry_list[0].fold, 2)

    def test_brute_force_reflection_symmetry_returns_list_of_two_reflection_planes(self):
        nuclei_array = MoleculeFactory.center_molecule(self.nuclei_array_h2o)
        rotation_symmetry_list = MoleculeFactory.brute_force_rotation_symmetry(nuclei_array)
        reflection_symmetry_list = MoleculeFactory.brute_force_reflection_symmetry(nuclei_array, rotation_symmetry_list)
        self.assertEqual(len(reflection_symmetry_list), 2)

    def test_check_linear_returns_false(self):
        nuclei_array = MoleculeFactory.center_molecule(self.nuclei_array_h2o)
        boolean = MoleculeFactory.check_linear(nuclei_array)
        self.assertEqual(boolean, False)

    def test_check_high_symmetry_returns_false(self):
        nuclei_array = MoleculeFactory.center_molecule(self.nuclei_array_h2o)
        symmetry_list = MoleculeFactory.brute_force_rotation_symmetry(nuclei_array)
        boolean = MoleculeFactory.check_high_symmetry(symmetry_list)
        self.assertEqual(boolean, False)


class TestSymmetryN2O(TestCase):

    def setUp(self):
        nitrogen_1 = MagicMock(element='NITROGEN', charge=7, mass=14, coordinates=(0.0000000000, 0.0000000000, -2.2684205883))
        nitrogen_2 = MagicMock(element='NITROGEN', charge=7, mass=14, coordinates=(0.0000000000, 0.0000000000, -0.1349300877))
        oxygen_1 = MagicMock(element='OXYGEN', charge=8, mass=16, coordinates=(0.0000000000, 0.0000000000, 2.1042369647))
        self.nuclei_array_n2o = [nitrogen_1, nitrogen_2, oxygen_1]

    def test_check_linear_returns_true(self):
        nuclei_array = MoleculeFactory.center_molecule(self.nuclei_array_n2o)
        boolean = MoleculeFactory.check_linear(nuclei_array)
        self.assertEqual(boolean, True)

    def test_check_inversion_symmetry_returns_true(self):
        nuclei_array = MoleculeFactory.center_molecule(self.nuclei_array_n2o)
        boolean = MoleculeFactory.check_inversion_symmetry(nuclei_array)
        self.assertEqual(boolean, False)

    def test_point_group_returns_c_inf_v_symmetry_for_nitrous_oxide(self):
        symmetry = MoleculeFactory.point_group(self.nuclei_array_n2o).symmetry_group
        testing.assert_equal(symmetry, 'C_{inf v}')


class TestSymmetryN2(TestCase):

    def setUp(self):
        nitrogen_1 = MagicMock(element='NITROGEN', charge=7, mass=14, coordinates=(0.0000000000, 0.0000000000, 1.0399092291))
        nitrogen_2 = MagicMock(element='NITROGEN', charge=7, mass=14, coordinates=(0.0000000000, 0.0000000000, -1.0399092291))
        self.nuclei_array_n2 = [nitrogen_1, nitrogen_2]

    def test_check_linear_returns_true(self):
        nuclei_array = MoleculeFactory.center_molecule(self.nuclei_array_n2)
        boolean = MoleculeFactory.check_linear(nuclei_array)
        self.assertEqual(boolean, True)

    def test_check_inversion_symmetry_returns_true(self):
        nuclei_array = MoleculeFactory.center_molecule(self.nuclei_array_n2)
        boolean = MoleculeFactory.check_inversion_symmetry(nuclei_array)
        self.assertEqual(boolean, True)

    def test_point_group_returns_d_inf_h_symmetry_for_nitrogen(self):
        symmetry = MoleculeFactory.point_group(self.nuclei_array_n2).symmetry_group
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

    def test_brute_force_rotation_symmetry_returns_list_of_four_axis_of_rotations_with_n_three(self):
        nuclei_array = MoleculeFactory.center_molecule(self.nuclei_array_ch4)
        symmetry_list = MoleculeFactory.brute_force_rotation_symmetry(nuclei_array)
        self.assertEqual([symmetry.fold for symmetry in symmetry_list].count(3), 4)

    def test_brute_force_rotation_symmetry_returns_list_of_four_axis_of_rotations_with_n_two(self):
        nuclei_array = MoleculeFactory.center_molecule(self.nuclei_array_ch4)
        symmetry_list = MoleculeFactory.brute_force_rotation_symmetry(nuclei_array)
        self.assertEqual([symmetry.fold for symmetry in symmetry_list].count(2), 3)

    def test_brute_force_reflection_symmetry_returns_a_list_of_six_reflection_planes(self):
        nuclei_array = MoleculeFactory.center_molecule(self.nuclei_array_ch4)
        rotation_symmetry_list = MoleculeFactory.brute_force_rotation_symmetry(nuclei_array)
        reflection_symmetry_list = MoleculeFactory.brute_force_reflection_symmetry(nuclei_array, rotation_symmetry_list)
        self.assertEqual(len(reflection_symmetry_list), 6)

    def test_check_linear_returns_false(self):
        nuclei_array = MoleculeFactory.center_molecule(self.nuclei_array_ch4)
        boolean = MoleculeFactory.check_linear(nuclei_array)
        self.assertEqual(boolean, False)

    def test_check_high_symmetry_returns_true(self):
        nuclei_array = MoleculeFactory.center_molecule(self.nuclei_array_ch4)
        symmetry_list = MoleculeFactory.brute_force_rotation_symmetry(nuclei_array)
        boolean = MoleculeFactory.check_high_symmetry(symmetry_list)
        self.assertEqual(boolean, True)

    def test_check_inversion_symmetry_returns_false(self):
        nuclei_array = MoleculeFactory.center_molecule(self.nuclei_array_ch4)
        boolean = MoleculeFactory.check_inversion_symmetry(nuclei_array)
        self.assertEqual(boolean, False)

    def test_point_group_returns_t_d_symmetry_for_methane(self):
        symmetry = MoleculeFactory.point_group(self.nuclei_array_ch4).symmetry_group
        testing.assert_equal(symmetry, 'T_{d}')
