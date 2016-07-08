from unittest import TestCase
from unittest.mock import MagicMock
from numpy import testing
from src.main.objects import MoleculeFactory


class TestSymmetryC2H4(TestCase):

    def setUp(self):
        carbon_1 = MagicMock(element='CARBON', charge=6, mass=12, coordinates=(0.0000000000, 1.2594652672, 0.0000000000))
        carbon_2 = MagicMock(element='CARBON', charge=6, mass=12, coordinates=(0.0000000000, -1.2594652672, 0.0000000000))
        hydrogen_1 = MagicMock(element='HYDROGEN', charge=1, mass=1, coordinates=(1.7400646600, 2.3216269636, 0.0000000000))
        hydrogen_2 = MagicMock(element='HYDROGEN', charge=1, mass=1, coordinates=(-1.7400646600, 2.3216269636, 0.0000000000))
        hydrogen_3 = MagicMock(element='HYDROGEN', charge=1, mass=1, coordinates=(1.7400646600, -2.3216269636, 0.0000000000))
        hydrogen_4 = MagicMock(element='HYDROGEN', charge=1, mass=1, coordinates=(-1.7400646600, -2.3216269636, 0.0000000000))
        self.nuclei_array_c2h4 = [carbon_1, carbon_2, hydrogen_1, hydrogen_2, hydrogen_3, hydrogen_4]
        self.molecule_factory = MoleculeFactory()

    def test_brute_force_rotation_symmetry_returns_list_of_three_axis_of_rotations(self):
        nuclei_array = self.molecule_factory.center_molecule(self.nuclei_array_c2h4)
        rotation, reflection, inversion = self.molecule_factory.brute_force_symmetry(nuclei_array)
        self.assertEqual(len(rotation), 3)

    def test_brute_force_reflection_symmetry_returns_list_of_three_reflection_planes(self):
        nuclei_array = self.molecule_factory.center_molecule(self.nuclei_array_c2h4)
        rotation, reflection, inversion = self.molecule_factory.brute_force_symmetry(nuclei_array)
        self.assertEqual(len(reflection), 3)

    def test_check_linear_returns_false(self):
        nuclei_array = self.molecule_factory.center_molecule(self.nuclei_array_c2h4)
        rotation, reflection, inversion = self.molecule_factory.brute_force_symmetry(nuclei_array)
        self.molecule_factory.standard_orientation(nuclei_array, rotation, reflection)
        boolean = self.molecule_factory.check_linear(nuclei_array)
        self.assertEqual(boolean, False)

    def test_check_high_symmetry_returns_false(self):
        nuclei_array = self.molecule_factory.center_molecule(self.nuclei_array_c2h4)
        rotation, reflection, inversion = self.molecule_factory.brute_force_symmetry(nuclei_array)
        self.molecule_factory.standard_orientation(nuclei_array, rotation, reflection)
        boolean = self.molecule_factory.check_high_symmetry(rotation)
        self.assertEqual(boolean, False)

    def test_check_n_two_fold_rotation_perpendicular_to_n_fold_returns_true(self):
        nuclei_array = self.molecule_factory.center_molecule(self.nuclei_array_c2h4)
        rotation, reflection, inversion = self.molecule_factory.brute_force_symmetry(nuclei_array)
        self.molecule_factory.standard_orientation(nuclei_array, rotation, reflection)
        boolean = self.molecule_factory.check_n_two_fold_perpendicular_to_n_fold(rotation)
        self.assertEqual(boolean, True)

    def test_check_sigma_h_returns_true(self):
        nuclei_array = self.molecule_factory.center_molecule(self.nuclei_array_c2h4)
        rotation, reflection, inversion = self.molecule_factory.brute_force_symmetry(nuclei_array)
        self.molecule_factory.standard_orientation(nuclei_array, rotation, reflection)
        boolean = self.molecule_factory.check_sigma_h(reflection)
        self.assertEqual(boolean, True)

    def test_point_group_returns_d_2h_symmetry_for_ethene(self):
        symmetry = self.molecule_factory.point_group(self.nuclei_array_c2h4).point_group.label
        testing.assert_equal(symmetry, 'D_{2h}')
