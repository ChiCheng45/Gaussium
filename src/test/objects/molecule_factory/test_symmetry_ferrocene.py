from unittest import TestCase
from unittest.mock import MagicMock
from numpy import testing
from src.main.objects import MoleculeFactory


class TestSymmetryFerrocene(TestCase):

    def setUp(self):
        carbon_1 = MagicMock(element='CARBON', charge=6, mass=12, coordinates=(0.99234, 0.72098, -1.73785))
        carbon_2 = MagicMock(element='CARBON', charge=6, mass=12, coordinates=(0.99234, -0.72098, -1.73785))
        carbon_3 = MagicMock(element='CARBON', charge=6, mass=12, coordinates=(-0.37904, -1.16657, -1.73785))
        carbon_4 = MagicMock(element='CARBON', charge=6, mass=12, coordinates=(-1.22660, 0.00000, -1.73785))
        carbon_5 = MagicMock(element='CARBON', charge=6, mass=12, coordinates=(-0.37904, 1.16657, -1.73785))
        hydrogen_1 = MagicMock(element='HYDROGEN', charge=1, mass=1, coordinates=(1.86677, 1.35629, -1.73016))
        hydrogen_2 = MagicMock(element='HYDROGEN', charge=1, mass=1, coordinates=(1.86677, -1.35629, -1.73016))
        hydrogen_3 = MagicMock(element='HYDROGEN', charge=1, mass=1, coordinates=(-0.71304, -2.19452, -1.73016))
        hydrogen_4 = MagicMock(element='HYDROGEN', charge=1, mass=1, coordinates=(-2.30745, 0.00000, -1.73016))
        hydrogen_5 = MagicMock(element='HYDROGEN', charge=1, mass=1, coordinates=(-0.71304, 2.19452, -1.73016))
        iron_1 = MagicMock(element='IRON', charge=26, mass=56, coordinates=(-0.00000, -0.00000, -0.00000))
        carbon_6 = MagicMock(element='CARBON', charge=6, mass=12, coordinates=(-0.99234, -0.72098, 1.73785))
        carbon_7 = MagicMock(element='CARBON', charge=6, mass=12, coordinates=(-0.99234, 0.72098, 1.73785))
        carbon_8 = MagicMock(element='CARBON', charge=6, mass=12, coordinates=(0.37904, 1.16657, 1.73785))
        carbon_9 = MagicMock(element='CARBON', charge=6, mass=12, coordinates=(1.22660, -0.00000, 1.73785))
        carbon_10 = MagicMock(element='CARBON', charge=6, mass=12, coordinates=(0.37904, -1.16657, 1.73785))
        hydrogen_6 = MagicMock(element='HYDROGEN', charge=1, mass=1, coordinates=(-1.86677, -1.35629, 1.73016))
        hydrogen_7 = MagicMock(element='HYDROGEN', charge=1, mass=1, coordinates=(-1.86677, 1.35629, 1.73016))
        hydrogen_8 = MagicMock(element='HYDROGEN', charge=1, mass=1, coordinates=(0.71304, 2.19452, 1.73016))
        hydrogen_9 = MagicMock(element='HYDROGEN', charge=1, mass=1, coordinates=(2.30745, -0.00000, 1.73016))
        hydrogen_10 = MagicMock(element='HYDROGEN', charge=1, mass=1, coordinates=(0.71304, -2.19452, 1.73016))
        self.nuclei_array_ferrocene = [carbon_1, carbon_2, carbon_3, carbon_4, carbon_5, carbon_6, carbon_7, carbon_8,
        carbon_9, carbon_10, hydrogen_1, hydrogen_2, hydrogen_3, hydrogen_4, hydrogen_5, hydrogen_6, hydrogen_7,
        hydrogen_8, hydrogen_9, hydrogen_10, iron_1]
        self.molecule_factory = MoleculeFactory()

    def test_point_group_returns_d_5d_symmetry_for_cubane(self):
        symmetry = self.molecule_factory.point_group(self.nuclei_array_ferrocene).point_group
        testing.assert_equal(symmetry, 'D_{5d}')
