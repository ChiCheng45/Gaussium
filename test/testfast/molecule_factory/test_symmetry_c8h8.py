from unittest import TestCase
from unittest.mock import MagicMock

from numpy import testing

from src.factory import MoleculeFactory
from src.factory import SymmetryFactory
from src.objects import InversionSymmetry


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
        self.molecule_factory = MoleculeFactory(symmetry=True)
        self.symmetry_factory = SymmetryFactory()
        self.inversion_symmetry = InversionSymmetry()

    def test_check_linear_returns_false(self):
        nuclei_array = self.molecule_factory.center_molecule(self.nuclei_array_c8h8)
        rotation, reflection, improper, inversion = self.symmetry_factory.brute_force_symmetry(nuclei_array)
        self.molecule_factory.standard_orientation(nuclei_array, rotation, reflection)
        boolean = self.molecule_factory.check_linear(nuclei_array)
        self.assertEqual(boolean, False)

    def test_check_high_symmetry_returns_true(self):
        nuclei_array = self.molecule_factory.center_molecule(self.nuclei_array_c8h8)
        rotation, reflection, improper, inversion = self.symmetry_factory.brute_force_symmetry(nuclei_array)
        self.molecule_factory.standard_orientation(nuclei_array, rotation, reflection)
        boolean = self.molecule_factory.check_high_symmetry(rotation)
        self.assertEqual(boolean, True)

    def test_point_group_returns_o_h_symmetry_for_cubane(self):
        symmetry = self.molecule_factory.create(self.nuclei_array_c8h8)[1].label
        testing.assert_equal(symmetry, 'O_{h}')
