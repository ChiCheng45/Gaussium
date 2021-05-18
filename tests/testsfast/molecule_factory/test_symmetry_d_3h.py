from unittest import TestCase
from unittest.mock import MagicMock
from numpy import testing
from gaussium.factory import MoleculeFactory
from gaussium.factory import SymmetryFactory


class TestSymmetryD3H(TestCase):

    def setUp(self):
        particle_1 = MagicMock(element='HYDROGEN', charge=1, mass=1, coordinates=(0.4391356726, 0.1106588251, -0.4635601962))
        particle_2 = MagicMock(element='HYDROGEN', charge=1, mass=1, coordinates=(-0.5185079933, 0.3850176090, 0.0537084789))
        particle_3 = MagicMock(element='HYDROGEN', charge=1, mass=1, coordinates=(0.0793723207, -0.4956764341, 0.4098517173))
        self.nuclei_array = [particle_1, particle_2, particle_3]
        self.molecule_factory = MoleculeFactory(symmetry=True)
        self.symmetry_factory = SymmetryFactory()

    def test_point_group_returns_d_3h_symmetry_for_system(self):
        symmetry = self.molecule_factory.create(self.nuclei_array)[1].label
        testing.assert_equal(symmetry, 'D_{3h}')
