from unittest import TestCase
from unittest.mock import MagicMock
from numpy import testing
from src.main.factory import MoleculeFactory
from src.main.factory import SymmetryFactory


class TestSymmetryN2O(TestCase):

    def setUp(self):
        nitrogen_1 = MagicMock(element='NITROGEN', charge=7, mass=14, coordinates=(0.0000000000, 0.0000000000, -2.2684205883))
        nitrogen_2 = MagicMock(element='NITROGEN', charge=7, mass=14, coordinates=(0.0000000000, 0.0000000000, -0.1349300877))
        oxygen_1 = MagicMock(element='OXYGEN', charge=8, mass=16, coordinates=(0.0000000000, 0.0000000000, 2.1042369647))
        self.nuclei_array_n2o = [nitrogen_1, nitrogen_2, oxygen_1]
        self.molecule_factory = MoleculeFactory(symmetry=True)
        self.symmetry_factory = SymmetryFactory()

    def test_check_linear_returns_true(self):
        nuclei_array = self.molecule_factory.center_molecule(self.nuclei_array_n2o)
        rotation, reflection, improper, inversion = self.symmetry_factory.brute_force_symmetry(nuclei_array)
        self.molecule_factory.standard_orientation(nuclei_array, rotation, reflection)
        boolean = self.molecule_factory.check_linear(nuclei_array)
        self.assertEqual(boolean, True)

    def test_point_group_returns_c_inf_v_symmetry_for_nitrous_oxide(self):
        symmetry = self.molecule_factory.create(self.nuclei_array_n2o).point_group.label
        testing.assert_equal(symmetry, 'C_{4v}')
