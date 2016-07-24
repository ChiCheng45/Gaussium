from unittest import TestCase
from unittest.mock import MagicMock
from numpy import testing
from src.main.objects import MoleculeFactory
from src.main.objects import SymmetryFactory


class TestSymmetryC2H6(TestCase):

    def setUp(self):
        carbon_1 = MagicMock(element='CARBON', charge=6, mass=12, coordinates=(0.7516, -0.0225, -0.0209))
        carbon_2 = MagicMock(element='CARBON', charge=6, mass=12, coordinates=(-0.7516, 0.0225, 0.0209))
        hydrogen_1 = MagicMock(element='HYDROGEN', charge=1, mass=1, coordinates=(1.1851, -0.0039, 0.9875))
        hydrogen_2 = MagicMock(element='HYDROGEN', charge=1, mass=1, coordinates=(1.1669, 0.8330, -0.5693))
        hydrogen_3 = MagicMock(element='HYDROGEN', charge=1, mass=1, coordinates=(1.1155, -0.9329, -0.5145))
        hydrogen_4 = MagicMock(element='HYDROGEN', charge=1, mass=1, coordinates=(-1.1669, -0.8334, 0.5687))
        hydrogen_5 = MagicMock(element='HYDROGEN', charge=1, mass=1, coordinates=(-1.1157, 0.9326, 0.5151))
        hydrogen_6 = MagicMock(element='HYDROGEN', charge=1, mass=1, coordinates=(-1.1850, 0.0044, -0.9875))
        self.nuclei_array_c2h6 = [carbon_1, carbon_2, hydrogen_1, hydrogen_2, hydrogen_3, hydrogen_4, hydrogen_5,
        hydrogen_6]
        self.molecule_factory = MoleculeFactory()
        self.symmetry_factory = SymmetryFactory()

    def test_check_linear_returns_false(self):
        nuclei_array = self.molecule_factory.center_molecule(self.nuclei_array_c2h6)
        rotation, reflection, improper, inversion = self.symmetry_factory.brute_force_symmetry(nuclei_array)
        self.molecule_factory.standard_orientation(nuclei_array, rotation, reflection)
        boolean = self.molecule_factory.check_linear(nuclei_array)
        self.assertEqual(boolean, False)

    def test_check_high_symmetry_returns_false(self):
        nuclei_array = self.molecule_factory.center_molecule(self.nuclei_array_c2h6)
        rotation, reflection, improper, inversion = self.symmetry_factory.brute_force_symmetry(nuclei_array)
        self.molecule_factory.standard_orientation(nuclei_array, rotation, reflection)
        boolean = self.molecule_factory.check_high_symmetry(rotation)
        self.assertEqual(boolean, False)

    def test_check_n_two_fold_rotation_perpendicular_to_n_fold_returns_true(self):
        nuclei_array = self.molecule_factory.center_molecule(self.nuclei_array_c2h6)
        rotation, reflection, improper, inversion = self.symmetry_factory.brute_force_symmetry(nuclei_array)
        self.molecule_factory.standard_orientation(nuclei_array, rotation, reflection)
        boolean = self.molecule_factory.check_n_two_fold_perpendicular_to_n_fold(rotation)
        self.assertEqual(boolean, True)

    def test_check_sigma_h_returns_false(self):
        nuclei_array = self.molecule_factory.center_molecule(self.nuclei_array_c2h6)
        rotation, reflection, improper, inversion = self.symmetry_factory.brute_force_symmetry(nuclei_array)
        self.molecule_factory.standard_orientation(nuclei_array, rotation, reflection)
        boolean = self.molecule_factory.check_sigma_h(reflection)
        self.assertEqual(boolean, False)

    def test_check_n_sigma_v_returns_true(self):
        nuclei_array = self.molecule_factory.center_molecule(self.nuclei_array_c2h6)
        rotation, reflection, improper, inversion = self.symmetry_factory.brute_force_symmetry(nuclei_array)
        self.molecule_factory.standard_orientation(nuclei_array, rotation, reflection)
        boolean = self.molecule_factory.check_n_sigma_v(3, reflection)
        self.assertEqual(boolean, True)

    def test_point_group_returns_d_3d_symmetry_for_cubane(self):
        symmetry = self.molecule_factory.point_group(self.nuclei_array_c2h6).point_group.label
        testing.assert_equal(symmetry, 'D_{3d}')
