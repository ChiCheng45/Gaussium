from unittest import TestCase
from unittest.mock import MagicMock

from numpy import testing

from src.factory import MoleculeFactory
from src.factory import SymmetryFactory


class TestSymmetryN2(TestCase):

    def setUp(self):
        nitrogen_1 = MagicMock(element='NITROGEN', charge=7, mass=14, coordinates=(0.0000000000, 0.0000000000, 1.0399092291))
        nitrogen_2 = MagicMock(element='NITROGEN', charge=7, mass=14, coordinates=(0.0000000000, 0.0000000000, -1.0399092291))
        self.nuclei_array_n2 = [nitrogen_1, nitrogen_2]
        self.molecule_factory = MoleculeFactory(symmetry=True)
        self.symmetry_factory = SymmetryFactory()

    def test_check_linear_returns_true(self):
        nuclei_array = self.molecule_factory.center_molecule(self.nuclei_array_n2)
        rotation, reflection, improper, inversion = self.symmetry_factory.brute_force_symmetry(nuclei_array)
        self.molecule_factory.standard_orientation(nuclei_array, rotation, reflection)
        boolean = self.molecule_factory.check_linear(nuclei_array)
        self.assertEqual(boolean, True)

    # D_{4h} because D_{inf h} is not useful for reducing integrals
    def test_point_group_returns_d_4_h_symmetry_for_nitrogen(self):
        symmetry = self.molecule_factory.create(self.nuclei_array_n2)[1].label
        testing.assert_equal(symmetry, 'D_{4h}')
