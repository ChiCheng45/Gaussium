from unittest import TestCase
from unittest.mock import MagicMock
from numpy import testing
from src.main.objects import MoleculeFactory


class TestSymmetryN2(TestCase):

    def setUp(self):
        nitrogen_1 = MagicMock(element='NITROGEN', charge=7, mass=14, coordinates=(0.0000000000, 0.0000000000, 1.0399092291))
        nitrogen_2 = MagicMock(element='NITROGEN', charge=7, mass=14, coordinates=(0.0000000000, 0.0000000000, -1.0399092291))
        self.nuclei_array_n2 = [nitrogen_1, nitrogen_2]
        self.molecule_factory = MoleculeFactory()

    def test_check_linear_returns_true(self):
        nuclei_array = self.molecule_factory.center_molecule(self.nuclei_array_n2)
        rotation, reflection = self.molecule_factory.brute_force_symmetry(nuclei_array)
        nuclei_array, rotation, reflection = self.molecule_factory.standard_orientation(nuclei_array, rotation, reflection)
        boolean = self.molecule_factory.check_linear(nuclei_array)
        self.assertEqual(boolean, True)

    def test_check_inversion_symmetry_returns_true(self):
        nuclei_array = self.molecule_factory.center_molecule(self.nuclei_array_n2)
        rotation, reflection = self.molecule_factory.brute_force_symmetry(nuclei_array)
        nuclei_array, rotation, reflection = self.molecule_factory.standard_orientation(nuclei_array, rotation, reflection)
        boolean = self.molecule_factory.check_inversion_symmetry(nuclei_array)
        self.assertEqual(boolean, True)

    def test_point_group_returns_d_inf_h_symmetry_for_nitrogen(self):
        symmetry = self.molecule_factory.point_group(self.nuclei_array_n2).point_group
        testing.assert_equal(symmetry, 'D_{inf h}')