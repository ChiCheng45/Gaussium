from unittest import TestCase
from unittest.mock import MagicMock
from numpy import testing
from src.main.objects import MoleculeFactory
from src.main.objects import SymmetryFactory


class TestSymmetryH2O(TestCase):

    def setUp(self):
        oxygen_1 = MagicMock(element='OXYGEN', charge=8, mass=16, coordinates=(0.0000000000, 0.0000000000, -0.1363928482))
        hydrogen_1 = MagicMock(element='HYDROGEN', charge=1, mass=1, coordinates=(0.0000000000, 1.4236595095, 0.9813433754))
        hydrogen_2 = MagicMock(element='HYDROGEN', charge=1, mass=1, coordinates=(0.0000000000, -1.4236595095, 0.9813433754))
        self.nuclei_array_h2o = [oxygen_1, hydrogen_1, hydrogen_2]
        self.molecule_factory = MoleculeFactory()
        self.symmetry_factory = SymmetryFactory()

    def test_brute_force_rotation_symmetry_returns_list_of_one_axis_of_rotations(self):
        nuclei_array = self.molecule_factory.center_molecule(self.nuclei_array_h2o)
        rotation, reflection, improper, inversion = self.symmetry_factory.brute_force_symmetry(nuclei_array)
        self.assertEqual(len(rotation), 1)

    def test_brute_force_rotation_symmetry_returns_axis_of_rotation_of_two_fold_symmetry(self):
        nuclei_array = self.molecule_factory.center_molecule(self.nuclei_array_h2o)
        rotation, reflection, improper, inversion = self.symmetry_factory.brute_force_symmetry(nuclei_array)
        self.assertEqual(rotation[0].fold, 2)

    def test_brute_force_reflection_symmetry_returns_list_of_two_reflection_planes(self):
        nuclei_array = self.molecule_factory.center_molecule(self.nuclei_array_h2o)
        rotation, reflection, improper, inversion = self.symmetry_factory.brute_force_symmetry(nuclei_array)
        self.assertEqual(len(reflection), 2)

    def test_check_linear_returns_false(self):
        nuclei_array = self.molecule_factory.center_molecule(self.nuclei_array_h2o)
        rotation, reflection, improper, inversion = self.symmetry_factory.brute_force_symmetry(nuclei_array)
        self.molecule_factory.standard_orientation(nuclei_array, rotation, reflection)
        boolean = self.molecule_factory.check_linear(nuclei_array)
        self.assertEqual(boolean, False)

    def test_check_high_symmetry_returns_false(self):
        nuclei_array = self.molecule_factory.center_molecule(self.nuclei_array_h2o)
        rotation, reflection, improper, inversion = self.symmetry_factory.brute_force_symmetry(nuclei_array)
        self.molecule_factory.standard_orientation(nuclei_array, rotation, reflection)
        boolean = self.molecule_factory.check_high_symmetry(rotation)
        self.assertEqual(boolean, False)

    def test_get_n_fold_returns_two(self):
        nuclei_array = self.molecule_factory.center_molecule(self.nuclei_array_h2o)
        rotation, reflection, improper, inversion = self.symmetry_factory.brute_force_symmetry(nuclei_array)
        self.molecule_factory.standard_orientation(nuclei_array, rotation, reflection)
        n = self.molecule_factory.return_principal_axis(rotation).fold
        self.assertEqual(n, 2)

    def test_check_n_two_fold_rotation_perpendicular_to_n_fold_returns_false(self):
        nuclei_array = self.molecule_factory.center_molecule(self.nuclei_array_h2o)
        rotation, reflection, improper, inversion = self.symmetry_factory.brute_force_symmetry(nuclei_array)
        self.molecule_factory.standard_orientation(nuclei_array, rotation, reflection)
        boolean = self.molecule_factory.check_n_two_fold_perpendicular_to_n_fold(rotation)
        self.assertEqual(boolean, False)

    def test_check_sigma_h_returns_false(self):
        nuclei_array = self.molecule_factory.center_molecule(self.nuclei_array_h2o)
        rotation, reflection, improper, inversion = self.symmetry_factory.brute_force_symmetry(nuclei_array)
        self.molecule_factory.standard_orientation(nuclei_array, rotation, reflection)
        boolean = self.molecule_factory.check_sigma_h(reflection)
        self.assertEqual(boolean, False)

    def test_check_n_sigma_v_returns_true(self):
        nuclei_array = self.molecule_factory.center_molecule(self.nuclei_array_h2o)
        rotation, reflection, improper, inversion = self.symmetry_factory.brute_force_symmetry(nuclei_array)
        self.molecule_factory.standard_orientation(nuclei_array, rotation, reflection)
        boolean = self.molecule_factory.check_n_sigma_v(2, reflection)
        self.assertEqual(boolean, True)

    def test_point_group_returns_c_2v_symmetry_for_water(self):
        symmetry = self.molecule_factory.point_group(self.nuclei_array_h2o).point_group.label
        testing.assert_equal(symmetry, 'C_{2v}')
