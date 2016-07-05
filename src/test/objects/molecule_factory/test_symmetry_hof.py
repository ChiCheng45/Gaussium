from unittest import TestCase
from unittest.mock import MagicMock
from numpy import testing
from src.main.objects import MoleculeFactory


class TestSymmetryHOF(TestCase):

    def setUp(self):
        oxygen_1 = MagicMock(element='OXYGEN', charge=8, mass=16, coordinates=(-1.4186923158, 0.1090030362, 0.0000000000))
        hydrogen_1 = MagicMock(element='HYDROGEN', charge=1, mass=1, coordinates=(-1.7313653816, -1.6895740638, 0.0000000000))
        fluorine_1 = MagicMock(element='FLUORINE', charge=9, mass=19, coordinates=(1.2899273141, 0.0031592817, 0.0000000000))
        self.nuclei_array_hof = [oxygen_1, hydrogen_1, fluorine_1]
        self.molecule_factory = MoleculeFactory()

    def test_brute_force_rotation_symmetry_returns_list_of_zero_axis_of_rotations(self):
        nuclei_array = self.molecule_factory.center_molecule(self.nuclei_array_hof)
        rotation, reflection, inversion = self.molecule_factory.brute_force_symmetry(nuclei_array)
        self.assertEqual(len(rotation), 0)

    def test_brute_force_reflection_symmetry_returns_list_of_one_reflection_planes(self):
        nuclei_array = self.molecule_factory.center_molecule(self.nuclei_array_hof)
        rotation, reflection, inversion = self.molecule_factory.brute_force_symmetry(nuclei_array)
        self.assertEqual(len(reflection), 1)

    def test_check_linear_returns_false(self):
        nuclei_array = self.molecule_factory.center_molecule(self.nuclei_array_hof)
        rotation, reflection, inversion = self.molecule_factory.brute_force_symmetry(nuclei_array)
        self.molecule_factory.standard_orientation(nuclei_array, rotation, reflection)
        boolean = self.molecule_factory.check_linear(nuclei_array)
        self.assertEqual(boolean, False)

    def test_check_high_symmetry_returns_false(self):
        nuclei_array = self.molecule_factory.center_molecule(self.nuclei_array_hof)
        rotation, reflection, inversion = self.molecule_factory.brute_force_symmetry(nuclei_array)
        self.molecule_factory.standard_orientation(nuclei_array, rotation, reflection)
        boolean = self.molecule_factory.check_high_symmetry(rotation)
        self.assertEqual(boolean, False)

    def test_check_sigma_h_returns_true(self):
        nuclei_array = self.molecule_factory.center_molecule(self.nuclei_array_hof)
        rotation, reflection, inversion = self.molecule_factory.brute_force_symmetry(nuclei_array)
        self.molecule_factory.standard_orientation(nuclei_array, rotation, reflection)
        boolean = self.molecule_factory.check_sigma_h(reflection)
        self.assertEqual(boolean, True)

    def test_point_group_returns_c_s_symmetry_for_hypofluorous_acid(self):
        symmetry = self.molecule_factory.point_group(self.nuclei_array_hof).point_group.label
        testing.assert_equal(symmetry, 'C_{s}')
