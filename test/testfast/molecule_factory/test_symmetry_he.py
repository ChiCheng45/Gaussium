from unittest import TestCase
from unittest.mock import MagicMock

from numpy import testing

from src.factory import MoleculeFactory


class TestSymmetryHe(TestCase):

    def setUp(self):
        helium_1 = MagicMock(element='HELIUM', charge=2, mass=4, coordinates=(-0.98781, 0.41551, 0.00000))
        self.nuclei_array_he = [helium_1]
        self.molecule_factory = MoleculeFactory(symmetry=True)

    def test_move_nuclei_to_the_origin(self):
        helium = self.molecule_factory.create(self.nuclei_array_he)[0][0]
        testing.assert_array_equal(helium.coordinates, (0.0, 0.0, 0.0))

    # O_{h} to speed up integral evaluation
    def test_point_group_returns_D_4h_symmetry_for_helium(self):
        symmetry = self.molecule_factory.create(self.nuclei_array_he)[1].label
        testing.assert_equal(symmetry, 'O_{h}')
