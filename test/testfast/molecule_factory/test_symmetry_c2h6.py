from unittest import TestCase
from unittest.mock import MagicMock
from numpy import testing
from gaussium.factory import MoleculeFactory
from gaussium.factory import SymmetryFactory


class TestSymmetryC2H6(TestCase):

    def setUp(self):
        carbon_1 = MagicMock(element='CARBON', charge=6, mass=12, coordinates=(-0.5297, 0.0100, -0.5606))
        carbon_2 = MagicMock(element='CARBON', charge=6, mass=12, coordinates=(0.5297, -0.0100, 0.5606))
        hydrogen_1 = MagicMock(element='HYDROGEN', charge=1, mass=1, coordinates=(-0.1452, 0.5058, -1.4448))
        hydrogen_2 = MagicMock(element='HYDROGEN', charge=1, mass=1, coordinates=(-1.4207, 0.5374, -0.2390))
        hydrogen_3 = MagicMock(element='HYDROGEN', charge=1, mass=1, coordinates=(-0.8163, -0.9984, -0.8374))
        hydrogen_4 = MagicMock(element='HYDROGEN', charge=1, mass=1, coordinates=(0.1452, -0.5058, 1.4448))
        hydrogen_5 = MagicMock(element='HYDROGEN', charge=1, mass=1, coordinates=(0.8163, 0.9984, 0.8374))
        hydrogen_6 = MagicMock(element='HYDROGEN', charge=1, mass=1, coordinates=(1.4207, -0.5374, 0.2390))

        self.nuclei_array_c2h6 = [
            carbon_1, carbon_2, hydrogen_1, hydrogen_2, hydrogen_3, hydrogen_4, hydrogen_5, hydrogen_6
        ]
        self.molecule_factory = MoleculeFactory(symmetry=True)
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
        symmetry = self.molecule_factory.create(self.nuclei_array_c2h6)[1].label
        testing.assert_equal(symmetry, 'D_{3d}')
